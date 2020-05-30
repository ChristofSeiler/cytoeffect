#' Plot Point Estimates for MCLE Fit
#'
#' @aliases plot.cytoeffect_poisson_mcle
#' @method plot cytoeffect_poisson_mcle
#'
#' @import ggplot2
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom magrittr %>% %<>%
#' @import dplyr
#' @import tibble
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param x Object of class \code{cytoeffect_poisson} computed
#'   using \code{\link{poisson_lognormal}}
#' @param type A string with the parameter to plot:
#'   \code{type = "beta"}, \code{type = "sigma"}, or \code{type = "Cor"}
#' @param selection A vector of strings with a selection of protein names to plot
#' @param ... Other parameters
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' # fit = cytoeffect::poisson_lognormal_mcle(...)
#' # plot(fit)
plot.cytoeffect_poisson_mcle = function(x, type = "beta",
                                        selection = x$protein_names, ...) {

  tb_args = x$tb_args
  protein_names = x$protein_names
  condition = x$condition
  group = x$group

  if(type == "beta") {

    set.seed(0xdada)
    coefs = lapply(tb_args$fit, function(x) x$beta) %>% bind_cols %>% t
    colnames(coefs) = protein_names
    tb_args %>%
      bind_cols(as_tibble(coefs)) %>%
      tidyr::pivot_longer(cols = protein_names, names_to = "marker") %>%
      ggplot(aes_string("value", "marker", "color" = condition)) +
      geom_jitter(height = 0.2, alpha = 0.5) +
      theme(legend.position = "bottom", legend.title = element_blank()) +
      theme(axis.title.y = element_blank()) +
      xlab("log expected count") +
      ggtitle("Fixed Effects"~beta)

  } else if (type == "sigma") {

    set.seed(0xdada)
    sds = lapply(tb_args$fit, function(x) sqrt(diag(x$Sigma))) %>% bind_cols %>% t
    colnames(sds) = protein_names
    tb_args %>%
      bind_cols(as_tibble(sds)) %>%
      tidyr::pivot_longer(cols = protein_names, names_to = "marker") %>%
      ggplot(aes_string("value", "marker", "color" = condition)) +
      geom_jitter(height = 0.2, alpha = 0.5) +
      theme(legend.position = "bottom", legend.title = element_blank()) +
      theme(axis.title.y = element_blank()) +
      xlab("sigma") +
      ggtitle("Marker Standard Deviation"~sigma)

  } else if (type == "Cor") {

    Sigma_list = lapply(tb_args$fit, function(x) x$Sigma)
    Sigma_array = simplify2array(Sigma_list)
    Sigma_mean = apply(Sigma_array, 1:2, mean)
    colnames(Sigma_mean) = rownames(Sigma_mean) = protein_names
    ggcorrplot::ggcorrplot(cov2cor(Sigma_mean), hc.order = TRUE, type = "lower",
                           outline.color = "lightgray",
                           colors = c("#6D9EC1", "white", "#E46726")) +
      ggtitle(paste0("Marker Correlations (cell)")) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

  } else {

    stop("Plotting for this type is not implemented.")

  }

}
