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
#' @param ncores Number of cores
#' @param show_donors Include donor random effect
#' @param show_markers Include markers
#' @param repel Repel marker names
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' # fit = cytoeffect::poisson_lognormal(...)
#' # cytoeffect::plot_distatis(fit)
plot_distatis = function(obj, ncores = parallel::detectCores(),
                         show_donors = TRUE, show_markers = TRUE, repel = TRUE) {

  if (class(obj) != "cytoeffect_poisson")
    stop("Not a cytoeffect_poisson object.")

  arrow_color = "darkgray"
  marker_color = "darkgray"
  segment_color = "black"
  marker_size = 5
  seed = 0xdada

  # sample all tables
  sample_info_k = c(obj$group,obj$condition,"k")
  n_chains = length(obj$fit_mcmc@stan_args)
  stan_args = obj$fit_mcmc@stan_args[[1]]
  n_total_draws = n_chains * (stan_args$iter - stan_args$warmup)
  expr_median = mclapply(
    1:n_total_draws,
    function(i) {
      set.seed(seed)
      posterior_predictive_log_lambda(obj, k = i, show_donors = show_donors) %>%
        group_by(.dots = sample_info_k) %>%
        summarize_at(obj$protein_names,median)
      },
    mc.cores = ncores
    )

  # prepare three way arrary for distatis
  dist_matrix = lapply(expr_median, function(x) {
    as.matrix(dist(x[,-seq(sample_info_k)]))
  })
  dist_matrix_arr = array(
    0, dim = c(dim(dist_matrix[[1]]), length(dist_matrix))
    )
  for(i in 1:dim(dist_matrix_arr)[3]) {
    dist_matrix_arr[,,i] = dist_matrix[[i]]
  }

  # run distatis
  fit_distatis = DistatisR::distatis(dist_matrix_arr, nfact2keep = 2)
  distatis_coords = fit_distatis$res4Splus[["PartialF"]]
  consensus_coords = fit_distatis$res4Splus[["F"]]

  # reshape distatis results
  distatis_coords_list = lapply(
    1:dim(distatis_coords)[3],
    function(i) bind_cols(expr_median[[1]],
                          as_tibble(distatis_coords[,,1]))
  )
  tb_distatis_coords = distatis_coords_list %>%
    bind_rows() %>%
    ungroup()
  tb_distatis_coords %<>% dplyr::rename(MDS1 = `Factor 1`,
                                        MDS2 = `Factor 2`)
  tb_consensus_coords = bind_cols(
    expr_median[[1]][sample_info_k],
    as_tibble(consensus_coords)) %>%
    ungroup()
  tb_consensus_coords %<>% dplyr::rename(MDS1 = `Factor 1`,
                                         MDS2 = `Factor 2`)

  # prepare plot distatis canvas
  ggmds = ggplot(tb_distatis_coords, aes_string(x = "MDS1", y = "MDS2",
                                                color = obj$condition))

  # prepare circle of correlation data
  protein_sd = apply(as.data.frame(tb_distatis_coords)[,obj$protein_names],2,sd)
  # only keep makers that have some variability
  protein_selection = obj$protein_names[protein_sd != 0]
  # correlations between variables and MDS axes
  expr_cor = cor(as.data.frame(tb_distatis_coords)[,protein_selection],
                 tb_distatis_coords[,c("MDS1","MDS2")]) %>% as_tibble
  # scaling factor (otherwise too crowded)
  expr_cor %<>% add_column(protein_selection)
  # add arrows coordinates
  expr_cor %<>% add_column(x0 = rep(0,nrow(expr_cor)))
  expr_cor %<>% add_column(y0 = rep(0,nrow(expr_cor)))
  scale_arrow = max(sqrt(tb_consensus_coords$MDS1^2 +
                           tb_consensus_coords$MDS2^2))
  cor_max = max(sqrt(expr_cor$MDS1^2+expr_cor$MDS2^2))
  expr_cor %<>% mutate(
    MDS1 = scale_arrow * MDS1/cor_max,
    MDS2 = scale_arrow * MDS2/cor_max
    )

  # add uncertainty countours
  ggmds = ggmds +
    xlab("Factor 1") +
    ylab("Factor 2") +
    scale_color_manual(values = c("#5DA5DA", "#FAA43A"),
                       name = obj$condition) +
    scale_fill_manual(values = c("#5DA5DA", "#FAA43A"),
                      name = obj$condition) +
    stat_density_2d(aes_string(fill = obj$condition),
                    geom = "polygon", alpha = 0.05)

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

  # add line segments connecting donor centers
  con_levels = levels(pull(tb_consensus_coords, obj$condition))
  segments =
    left_join(
      tb_consensus_coords[tb_consensus_coords[,obj$condition] == con_levels[1],],
      tb_consensus_coords[tb_consensus_coords[,obj$condition] == con_levels[2],],
      by = obj$group
    )
  segments %<>% dplyr::select(
    obj$group, MDS1.x, MDS2.x, MDS1.y, MDS2.y
  )
  ggmds = ggmds + geom_segment(
    aes(x = MDS1.x, xend = segments$MDS1.y, y = MDS2.x, yend = MDS2.y),
    colour = segment_color, alpha = 1.0,
    data = segments)

  # add donors centers
  if(show_donors) {
    ggmds = ggmds +
      geom_point(
        data = tb_consensus_coords,
        aes_string(shape = obj$group,
                   x = tb_consensus_coords$MDS1,
                   y = tb_consensus_coords$MDS2),
        size = 4,
        color = "black"
        ) +
      scale_shape_manual(
        values =
          64 + # from shape table (so that it starts at A)
          pull(tb_consensus_coords, obj$group) %>%
          unique %>%
          length %>%
          seq(1, .)
      ) +
      theme(legend.position = "bottom")
  }

  # add title
  ggmds + ggtitle("Posterior DiSTATIS of Latent Variable"~lambda)

}
