#' Plot traceplots.
#'
#' @import rstan
#' @import ggplot2
#' @import reshape2
#' @export
#'
traceplot = function(obj) {

  if (class(obj) != "cytoeffect")
    stop("Not a cytoeffect object.")

  fit_mcmc = obj$fit_mcmc
  protein_names = obj$protein_names
  warmup = fit_mcmc@stan_args[[1]]$warmup

  beta_all = rstan::extract(fit_mcmc,
                            pars = c("beta[2]","beta[3]","beta[10]"),
                            permuted = FALSE,
                            inc_warmup = TRUE)
  beta_warmup = reshape2::melt(beta_all,
                               varnames = c("iteration","chain","parameter"))
  ggplot(beta_warmup, aes(x = iteration, y = value, color = parameter)) +
    annotate("rect",
             xmin = 0, xmax = warmup,
             ymin = -Inf, ymax = Inf,
             alpha = 0.2,
             color = "gray") +
    geom_line() +
    facet_wrap(~ chain)

}
