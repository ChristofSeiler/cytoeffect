#' Plot Multivariate Posterior Summaries for Poisson Log-Normal Mixed Model
#'
#' @import rstan
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import parallel
#' @import MASS
#' @export
#'
#' @param obj Object of class \code{cytoeffect_poisson} computed
#'   using \code{\link{poisson_lognormal}}
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
posterior_predictive_log_lambda = function(obj, k = 1, show_donors = TRUE) {

  if (class(obj) != "cytoeffect_poisson")
    stop("Not a cytoeffect_poisson object.")

  stan_pars = rstan::extract(obj$fit_mcmc,
                             pars = c("beta",
                                      "sigma","sigma_term",
                                      "Q","Q_term",
                                      "b_donor",
                                      "theta"))
  condition_index = seq(obj$conditions)
  donor = obj$df_samples_subset %>%
    pull(obj$group) %>%
    as.factor %>%
    levels
  donor_index = seq(donor)
  # kth posterior draw
  sample_condition_donor = function(k, tb_info) {
    # fixed effects
    beta = stan_pars$beta[k,,tb_info$term_index]
    # cell random effect
    if(tb_info$term_index == 1) {
      sigma = stan_pars$sigma[k,]
      Q = stan_pars$Q[k,,]
    } else {
      sigma = stan_pars$sigma_term[k,]
      Q = stan_pars$Q_term[k,,]
    }
    Cov = Q %*% diag(sigma^2) %*% t(Q)
    # combine
    if(show_donors) {
      # donor random effect
      b_donor = stan_pars$b_donor[k,tb_info$donor_index,]
      lambda = mvrnorm(n = tb_info$n, beta + b_donor, Cov)
    } else {
      lambda = mvrnorm(n = tb_info$n, beta, Cov)
    }
    # account for zero inflation
    theta = stan_pars$theta[k,,]
    zeros = matrix(rbinom(tb_info$n*length(obj$protein_names),
                          size = 1, # number of trials is 1 for Bernoulli
                          prob = theta # mixture proportion
                          ),
                   nrow = length(obj$protein_names),
                   ncol = tb_info$n) %>% t
    #tibble(sim = apply(zeros, 2, mean), theta, diff = sim-theta)
    lambda[which(zeros == 1, arr.ind = TRUE)] = 0
    lambda %<>% as.tibble
    names(lambda) = obj$protein_names
    lambda %<>% add_column(term  = tb_info$term)
    lambda %<>% add_column(donor  = tb_info$donor)
    lambda %<>% add_column(k  = k)
    lambda
  }
  # count number of cells per term and donor
  subgroups = obj$df_samples_subset %>%
    group_by_(term = obj$condition, donor = obj$group) %>%
    tally %>%
    ungroup
  subgroups %<>% mutate(term_index = subgroups %>%
                          pull(term) %>%
                          as.integer)
  subgroups %<>% mutate(donor_index = subgroups %>%
                          pull(donor) %>%
                          as.factor %>%
                          as.integer)
  # sample one table
  lapply(seq(nrow(subgroups)), function(i) {
    sample_condition_donor(k = k, tb_info = subgroups[i,])
    }) %>% bind_rows()

}
