#' Plot Multivariate Posterior Summaries for Poisson Log-Normal Mixed Model
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import dplyr
#' @import parallel
#' @importFrom MASS mvrnorm
#' @importFrom rstan extract
#' @export
#'
#' @param obj Object of class \code{cytoeffect_poisson} computed
#'   using \code{\link{poisson_lognormal}}
#' @param k Draw from HMC chain
#' @param show_donors Include donor random effect
#' @return \code{\link[tibble]{tibble}} object
#'
posterior_predictive_log_lambda = function(obj, k = 1, show_donors = TRUE) {

  if(!is(obj, "cytoeffect_poisson"))
    stop("Not a cytoeffect_poisson object.")

  stan_pars = rstan::extract(obj$fit_mcmc,
                             pars = c("beta",
                                      "sigma",
                                      "Cor",
                                      "b_donor",
                                      "theta"))

  condition = obj$df_samples_subset %>%
    pull(obj$condition) %>%
    as.factor %>%
    levels
  condition_index = seq(obj$conditions)

  group = obj$df_samples_subset %>%
    pull(obj$group) %>%
    as.factor %>%
    levels
  group_index = seq(group)

  # kth posterior draw
  sample_condition_donor = function(k, tb_info) {
    # fixed effects
    beta = stan_pars$beta[k,,tb_info$cond_index]
    # cell random effect
    sigma = stan_pars$sigma[k,]
    Cor = stan_pars$Cor[k,,]
    Cov = diag(sigma) %*% Cor %*% diag(sigma)
    # combine
    if(show_donors) {
      # donor random effect
      b_donor = stan_pars$b_donor[k,tb_info$group_index,]
      lambda = MASS::mvrnorm(n = tb_info$n, beta + b_donor, Cov)
    } else {
      lambda = MASS::mvrnorm(n = tb_info$n, beta, Cov)
    }
    # make sure it has the right dimension
    lambda %<>% matrix(nrow = tb_info$n, ncol = length(obj$protein_names))
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
    colnames(lambda) = obj$protein_names
    lambda %<>% as_tibble
    lambda %<>% add_column(!!(obj$condition) := pull(tb_info, obj$condition))
    lambda %<>% add_column(!!(obj$group) := pull(tb_info, obj$group))
    lambda %<>% add_column(k  = k)
    lambda
  }
  # count number of cells per term and donor
  subgroups = obj$df_samples_subset %>%
    group_by_at(c(obj$condition, obj$group)) %>%
    tally %>%
    ungroup
  subgroups %<>% mutate(cond_index = subgroups %>%
                          pull(obj$condition) %>%
                          as.integer)
  subgroups %<>% mutate(group_index = subgroups %>%
                          pull(obj$group) %>%
                          as.factor %>%
                          as.integer)
  # sample one table
  lapply(seq(nrow(subgroups)), function(i) {
    sample_condition_donor(k = k, tb_info = subgroups[i,])
    }) %>% bind_rows()

}
