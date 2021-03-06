#' DiSTATIS plot of posterior or bootstrap samples
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import dplyr
#' @import parallel
#' @import ggrepel
#' @importFrom DistatisR distatis
#' @importFrom MASS mvrnorm
#' @importFrom tidyr unnest
#' @export
#'
#' @param obj Object of class \code{cytoeffect_poisson} or \code{cytoeffect_poisson_mcle}
#'   computed by \code{\link{poisson_lognormal}} or \code{\link{poisson_lognormal_mcle}}
#' @param ndraws Number of posterior or bootstrap samples
#' @param ncores Number of cores
#' @param show_donors Include donor random effect
#' @param show_markers Include markers
#' @param repel Repel marker names
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' df = simulate_data()
#' str(df)
#' fit = poisson_lognormal(df,
#'                         protein_names = names(df)[3:ncol(df)],
#'                         condition = "condition",
#'                         group = "donor",
#'                         r_donor = 2,
#'                         warmup = 200, iter = 325, adapt_delta = 0.95,
#'                         num_chains = 1)
#' plot_distatis(fit, ndraws = 125)
#' fit = poisson_lognormal_mcle(df,
#'                              protein_names = names(df)[3:ncol(df)],
#'                              condition = "condition",
#'                              group = "donor",
#'                              ncores = 1)
#' plot_distatis(fit, ndraws = 125)
#' }
#'
plot_distatis = function(obj, ndraws = 1000, ncores = 1,
                         show_donors = TRUE, show_markers = TRUE, repel = TRUE) {

  if(!is(obj, "cytoeffect_poisson") & !is(obj, "cytoeffect_poisson_mcle"))
    stop("Not a cytoeffect_poisson or cytoeffect_poisson_mcle object.")

  arrow_color = "darkgray"
  marker_color = "darkgray"
  segment_color = "black"
  marker_size = 5
  options(mc.cores = ncores)

  # sample all tables
  if(is(obj, "cytoeffect_poisson")) {
    # posterior samples
    sample_info_k = c(obj$group,obj$condition,"k")
    n_chains = length(obj$fit_mcmc@stan_args)
    stan_args = obj$fit_mcmc@stan_args[[1]]
    nsamples = n_chains * (stan_args$iter - stan_args$warmup)
    ndraws = min(nsamples, ndraws)
    expr_median = mclapply(
      seq_len(ndraws),
      function(i) {
        posterior_predictive_log_lambda(obj, k = i, show_donors = show_donors) %>%
          group_by(.dots = sample_info_k) %>%
          summarize_at(obj$protein_names,median)
        }
      )
  } else {
    # parametric bootstrap samples
    boot = function(tb_args) {
      boot_list = lapply(1:nrow(tb_args), function(i) {
        Y_donor = MASS::mvrnorm(n = tb_args$n[i],
                                mu = tb_args$fit[[i]]$beta,
                                Sigma = tb_args$fit[[i]]$Sigma)
        colnames(Y_donor) = obj$protein_names
        Y_donor %>% as_tibble %>% summarize_all(median)
      })
      tb_args %>%
        dplyr::select(-Y, -fit, -n) %>%
        add_column(b = boot_list) %>%
        tidyr::unnest(b) %>%
        ungroup()
    }
    expr_median = mclapply(1:ndraws, function(i) boot(obj$tb_args))
  }

  # prepare three way arrary for distatis
  dist_matrix = lapply(expr_median, function(x) {
    x %>%
      ungroup() %>%
      select_at(obj$protein_names) %>%
      dist %>%
      as.matrix
  })
  dist_matrix_arr = simplify2array(dist_matrix)

  # run distatis
  fit_distatis = DistatisR::distatis(dist_matrix_arr, nfact2keep = 2)
  distatis_coords = fit_distatis$res4Splus[["PartialF"]]
  consensus_coords = fit_distatis$res4Splus[["F"]]

  # reshape distatis results
  distatis_coords_list = lapply(
    1:dim(distatis_coords)[3],
    function(i) bind_cols(expr_median[[i]],
                          as_tibble(distatis_coords[,,i]))
  )
  tb_distatis_coords = distatis_coords_list %>%
    bind_rows() %>%
    ungroup()
  tb_distatis_coords %<>% dplyr::rename(MDS1 = `Factor 1`,
                                        MDS2 = `Factor 2`)
  tb_consensus_coords = bind_cols(
    expr_median[[1]] %>% dplyr::select_at(vars(obj$group, obj$condition)),
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
                                      alpha = 1.0)
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
    aes(x = MDS1.x, xend = MDS1.y, y = MDS2.x, yend = MDS2.y),
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
  if(is(obj, "cytoeffect_poisson")) {
    ggmds + ggtitle("Posterior DiSTATIS of Latent Variable"~lambda)
  } else {
    ggmds + ggtitle("Parametric Bootstrap DiSTATIS of Latent Variable"~lambda)
  }

}
