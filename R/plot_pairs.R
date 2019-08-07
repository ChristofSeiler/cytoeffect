#' Posterior multivariate pairs plot for Poisson Log-Normal Mixed Model
#'
#' @import rstan
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import cowplot
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
#' # fit = cytoeffect::poisson_lognormal(...)
#' # cytoeffect::plot_pairs(fit,
#'                          marker1 = "pSTAT1",
#'                          marker2 = "pSTAT3",
#'                          marker3 = "pSTAT5")
plot_pairs = function(obj, marker1, marker2, marker3) {

  if (class(obj) != "cytoeffect_poisson")
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
  plot_diag = function(marker) {
    ggplot(tb_log_count, aes_string(marker, fill = obj$condition)) +
      geom_histogram(bins = 40, position = "identity", alpha = 0.5) +
      scale_fill_few()
  }
  plot_off_diag = function(marker1, marker2) {
    ggplot(tb_log_count, aes_string(marker1, marker2, color = obj$condition)) +
      geom_density2d() +
      scale_color_few()
  }
  ppair = plot_grid(
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
  plot_grid(ppair,
            get_legend(plot_diag(marker1) + theme(legend.position = "bottom")),
            ncol = 1,
            rel_heights = c(1, .1))

}
