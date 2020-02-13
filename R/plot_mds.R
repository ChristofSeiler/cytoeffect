#' Plot Multivariate Posterior Summaries for Poisson Log-Normal Mixed Model
#'
#' @import rstan
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import parallel
#' @import MASS
#' @import ggrepel
#' @export
#'
#' @param obj Object of class \code{cytoeffect_poisson} computed
#'   using \code{\link{poisson_lognormal}}
#' @param asp Set \code{asp = FALSE} to avoid scaling aspect ratio by eigenvalues
#' @param ncores Number of cores
#' @param thinning Number of posterior draws after thinning
#' @param show_donors Include donor random effect
#' @param show_markers Include markers
#' @param cor_scaling_factor Scaling factor for correlation arrows
#' @param repel Repel marker names
#' @param scale_x Scale x axis
#' @param scale_y Scale y axis
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' # fit = cytoeffect::poisson_lognormal(...)
#' # cytoeffect::plot_mds(fit, asp = FALSE)
plot_mds = function(obj, asp = TRUE, ncores = parallel::detectCores(), thinning = 100,
                    show_donors = TRUE, show_markers = TRUE, cor_scaling_factor = 1,
                    repel = TRUE, scale_x = 1, scale_y = 1) {

  if (class(obj) != "cytoeffect_poisson")
    stop("Not a cytoeffect_poisson object.")

  arrow_color = "darkgray"
  marker_color = "black"
  segment_color = "black"
  marker_size = 5
  seed = 0xdada

  # sample all tables
  sample_info_k = c(obj$group,obj$condition,"k")

  n_chains = length(obj$fit_mcmc@stan_args)
  stan_args = obj$fit_mcmc@stan_args[[1]]
  n_total_draws = n_chains * (stan_args$iter - stan_args$warmup)
  expr_median = mclapply(
    round(seq(1, n_total_draws, length.out = thinning)),
    function(i) {
      set.seed(seed)
      posterior_predictive_log_lambda(obj, k = i, show_donors = show_donors) %>%
        group_by(.dots = sample_info_k) %>%
        summarize_at(obj$protein_names,median)
      },
    mc.cores = ncores
    ) %>% bind_rows
  # classical MDS on all posterior draws
  dist_matrix = dist(expr_median[,-seq(sample_info_k)])
  set.seed(seed)
  mds_res = cmdscale(dist_matrix,eig = TRUE, k = 2) # k is the number of dim
  explained_var = (100*mds_res$eig[1:2]/sum(mds_res$eig)) %>% round(digits = 1)
  expr_median %<>% bind_cols(tibble(MDS1 = mds_res$points[,1],
                                    MDS2 = mds_res$points[,2]))
  # plot MDS
  ggmds = ggplot(expr_median, aes_string(x = "MDS1", y = "MDS2",
                                         color = obj$condition))

  # make circle of correlation plot
  protein_sd = apply(as.data.frame(expr_median)[,obj$protein_names],2,sd)
  # only keep makers that have some variability
  protein_selection = obj$protein_names[protein_sd != 0]
  # correlations between variables and MDS axes
  expr_cor = cor(as.data.frame(expr_median)[,protein_selection],
                 expr_median[,c("MDS1","MDS2")]) %>% as_tibble
  # scaling factor (otherwise too crowded)
  expr_cor = expr_cor * cor_scaling_factor
  expr_cor %<>% add_column(protein_selection)
  # add arrows coordinates
  expr_cor %<>% add_column(x0 = rep(0,nrow(expr_cor)))
  expr_cor %<>% add_column(y0 = rep(0,nrow(expr_cor)))

  # add correlation arrows
  if(show_markers) {
    ggmds = ggmds + annotate("segment",
                             x = expr_cor$x0, xend = expr_cor$MDS1,
                             y = expr_cor$y0, yend = expr_cor$MDS2,
                             colour = arrow_color,
                             alpha = 1.0,
                             arrow = arrow(type = "open", length = unit(0.03, "npc")))

    ## add marker names labels
    if(repel) {
      ggmds = ggmds + geom_text_repel(data = expr_cor,
                                      aes(x = MDS1, y = MDS2,
                                          label = protein_selection),
                                      color = marker_color,
                                      alpha = 1.0, seed = seed)
    } else {
      ggmds = ggmds + geom_text(data = expr_cor,
                                aes(x = MDS1, y = MDS2,
                                    label = protein_selection),
                                color = marker_color,
                                alpha = 1.0)
    }
  }

  # add uncertainty countours
  ggmds = ggmds +
    xlab(paste0("MDS1 (",explained_var[1],"%)")) +
    ylab(paste0("MDS2 (",explained_var[2],"%)")) +
    scale_color_manual(values = c("#5DA5DA", "#FAA43A"),
                       name = obj$condition) +
    scale_fill_manual(values = c("#5DA5DA", "#FAA43A"),
                      name = obj$condition) +
    stat_density_2d(aes_string(fill = obj$condition),
                    geom = "polygon", alpha = 0.05)

  # add median donors
  if(show_donors) {
    expr_median_donor = expr_median %>%
      group_by(.dots = c(obj$group, obj$condition)) %>%
      summarize_at(c("MDS1","MDS2"), median)
    expr_median_donor %<>% add_column(
      color = sapply(pull(expr_median_donor, obj$condition),
                     function(x) if(x == pull(expr_median_donor, obj$condition)[1])
                       "#5DA5DA" else "#FAA43A"))
    ggmds = ggmds +
      geom_point(
        data = expr_median_donor,
        aes_string(shape = obj$group,
                   x = expr_median_donor$MDS1,
                   y = expr_median_donor$MDS2),
        size = 3,
        ) +
      scale_shape_manual(
        values =
          64 + # from shape table (so that it starts at A)
          pull(expr_median_donor, obj$group) %>%
            unique %>%
            length %>%
            seq(1, .)
        ) +
      theme(legend.position = "bottom")
  }

  # add line segments connecting donor centers
  con_levels = levels(pull(expr_median_donor, obj$condition))
  segments =
    left_join(
      expr_median_donor[expr_median_donor[,obj$condition] == con_levels[1],],
      expr_median_donor[expr_median_donor[,obj$condition] == con_levels[2],],
      by = obj$group
    )
  segments %<>% dplyr::select(
    donor, MDS1.x, MDS2.x, MDS1.y, MDS2.y
  )
  ggmds = ggmds + geom_segment(
    aes(x = MDS1.x, xend = segments$MDS1.y, y = MDS2.x, yend = MDS2.y),
    colour = segment_color, alpha = 1.0,
    data = segments)

  # add title
  ggmds = ggmds +
    ggtitle("Posterior MDS of Latent Variable"~lambda~"(Aspect Ratio Unscaled)")

  ggmds = ggmds +
    scale_x_continuous(limits = range(expr_median$MDS1)*scale_x) +
    scale_y_continuous(limits = range(expr_median$MDS2)*scale_y)

  if(asp) {
    # change aspect ratio according to explained variance
    ggmds = ggmds +
      coord_fixed(ratio = explained_var[2] / explained_var[1]) +
      ggtitle("Posterior MDS of Latent Variable"~lambda)
  }
  ggmds

}
