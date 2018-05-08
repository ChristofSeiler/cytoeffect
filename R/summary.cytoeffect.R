#' Posterior summary.
#'
#' @import rstan
#' @import ggplot2
#' @import tibble
#' @export
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
