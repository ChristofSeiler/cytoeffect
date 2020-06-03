#' Posterior multivariate pairs plot for Poisson Log-Normal Mixed Model
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom cowplot plot_grid get_legend
#' @importFrom ggthemes scale_colour_few
#' @importFrom rstan extract
#' @export
#'
#' @param obj Object of class \code{cytoeffect_poisson} computed
#'   using \code{\link{poisson_lognormal}}
#' @param marker1 Name of first marker
#' @param marker2 Name of second marker
#' @param marker3 Name of third marker
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' df = simulate_data(n_cells = 10)
#' str(df)
#' fit = poisson_lognormal(df,
#'                         protein_names = names(df)[3:ncol(df)],
#'                         condition = "condition",
#'                         group = "donor",
#'                         r_donor = 2,
#'                         warmup = 200, iter = 325, adapt_delta = 0.95,
#'                         num_chains = 1)
#' plot_pairs(fit, "m01", "m02", "m03")
#'
plot_pairs = function(obj, marker1, marker2, marker3) {

  if (!is(obj, "cytoeffect_poisson"))
    stop("Not a cytoeffect_poisson object.")

  pSTAT1_index = which(obj$protein_names == marker1)
  pSTAT3_index = which(obj$protein_names == marker2)
  pSTAT5_index = which(obj$protein_names == marker3)
  post_beta = rstan::extract(obj$fit_mcmc, pars = "beta")[[1]]
  tb_log_count = bind_rows(
    tibble(
      term = levels(pull(obj$df_samples_subset, obj$condition))[1],
      pSTAT1 = post_beta[,pSTAT1_index,1],
      pSTAT3 = post_beta[,pSTAT3_index,1],
      pSTAT5 = post_beta[,pSTAT5_index,1]
    ),
    tibble(
      term = levels(pull(obj$df_samples_subset, obj$condition))[2],
      pSTAT1 = post_beta[,pSTAT1_index,2],
      pSTAT3 = post_beta[,pSTAT3_index,2],
      pSTAT5 = post_beta[,pSTAT5_index,2]
    )
  )
  names(tb_log_count) = c(obj$condition, marker1, marker2, marker3)
  plot_diag = function(marker) {
    ggplot(tb_log_count, aes_string(marker, fill = obj$condition)) +
      geom_histogram(bins = 40, position = "identity", alpha = 0.5) +
      ggthemes::scale_fill_few()
  }
  plot_off_diag = function(marker1, marker2) {
    ggplot(tb_log_count, aes_string(marker1, marker2, color = obj$condition)) +
      geom_density2d() +
      ggthemes::scale_color_few()
  }
  ppair = cowplot::plot_grid(
    plot_diag(marker1) + theme(legend.position = "none"),
    NULL,
    NULL,
    plot_off_diag(marker1,marker2) + theme(legend.position = "none"),
    plot_diag(marker2) + theme(legend.position = "none"),
    NULL,
    plot_off_diag(marker1,marker3) + theme(legend.position = "none"),
    plot_off_diag(marker2,marker3) + theme(legend.position = "none"),
    plot_diag(marker3) + theme(legend.position = "none"),
    ncol = 3
  )
  cowplot::plot_grid(ppair,
                     cowplot::get_legend(plot_diag(marker1) +
                                           theme(legend.position = "bottom")),
                     ncol = 1,
                     rel_heights = c(1, .1))

}
