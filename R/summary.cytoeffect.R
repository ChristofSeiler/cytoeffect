#' Posterior Summary for Poisson Log-Normal Mixed Model
#'
#' @import rstan
#' @import ggplot2
#' @import tibble
#' @export
#'
#' @param obj Object of class \code{cytoeffect_poisson} computed using \code{\link{poisson_lognormal}}
#' @param type A string with the variable name to plot: \code{par = "beta"}, \code{par = "sigma"}, or \code{par = "Cor"}
#' @return \code{\link[tibble]{tibble}} object
#'
summary.cytoeffect = function(obj, par) {

  if (class(obj) != "cytoeffect")
    stop("Not a cytoeffect object.")

  fit_mcmc = obj$fit_mcmc
  protein_names = obj$protein_names

  param = rstan::extract(fit_mcmc,pars = par)[[1]]
  param = param[,-1] # remove intercept
  tibble(protein_name = protein_names,
         median = apply(param,2,median),
         low = apply(param,2,function(x) quantile(x,probs = 0.025)),
         high = apply(param,2,function(x) quantile(x,probs = 0.975)))

}
