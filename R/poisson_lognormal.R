#' HMC Sampling for Poisson Log-Normal Mixed Model
#'
#' \code{poisson_lognormal} uses Hamiltonion Monte Carlo to sample form an extend
#' Poisson log-normal mixed model. Each cell and protein marker has its own rate
#' parameter following a linear model.
#'
#' @import rstan
#' @import reshape2
#' @import dplyr
#' @import batchtools
#' @export
#'
#' @param df_samples_subset Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param protein_names A vector of column names of protein to use in the analysis
#' @param condition The column name of the condition variable
#' @param group The column name of the group variable
#' @param r_donor Rank of the donor random effect covariance matrix
#' @param iter Number of iteration per chain for the HMC sampler
#' @param warmup Number of warm up steps per chain for the HMC sampler
#' @param num_chains Number of HMC chains to run in parallel
#' @param adapt_delta Parameter to control step size of numerical solver
#' @param seed Set seed for HMC sampler
#'                    (higher value means smaller step size, max is 1)
#' @return A list of class \code{cytoeffect_poisson} containing
#'   \item{fit_mcmc}{\code{\link[rstan]{rstan}} object}
#'   \item{protein_names}{input protein names}
#'   \item{condition}{input condition variable}
#'   \item{group}{input group names}
#'   \item{df_samples_subset}{input df_samples_subset table}
#'
poisson_lognormal = function(df_samples_subset,
                             protein_names,
                             condition,
                             group,
                             r_donor,
                             eta = 1,
                             iter = 325,
                             warmup = 200,
                             num_chains = 4,
                             adapt_delta = 0.8,
                             seed = 1) {

  # some checks
  if(sum(names(df_samples_subset) == condition) == 0)
    stop("condition column missing")
  if(sum(names(df_samples_subset) == group) == 0)
    stop("group column missing")
  if(nrow(df_samples_subset) == 0)
    stop("no observations")
  if(df_samples_subset %>% pull(condition) %>% nlevels != 2)
    stop("condition variables should have two levels")

  # prepare input data
  Y = df_samples_subset %>%
    ungroup() %>%
    dplyr::select(protein_names) %>%
    as.matrix()
  term = df_samples_subset %>%
    pull(condition) %>%
    as.factor() %>%
    as.integer()
  p = length(table(term))
  n = nrow(Y)
  d = ncol(Y)
  donor = df_samples_subset %>%
    pull(group) %>%
    as.factor() %>%
    as.integer()
  k = length(table(donor))
  stan_data = list(Y = t(Y), n = n, d = d, p = p,
                   k = k, donor = donor, term = term,
                   r_donor = r_donor,
                   eta = eta)

  # prepare starting point for sampler

  # log mean per condition (same as Poisson regression)
  beta = df_samples_subset %>%
    group_by_at(condition) %>%
    summarise_at(protein_names, function(x) log(mean(x))) %>%
    dplyr::select(protein_names) %>%
    t
  # set according to prior
  is.neginf = function(x) x == -Inf
  beta[is.neginf(beta)] = -7 # beta[j] ~ normal(0, 7);

  # covariance matrix across cells per level of condition
  tfm = function(x) asinh(x/5)
  # transform raw counts
  Y_tfm = tfm(Y)
  # sample standard deviation
  Y_cov = cov(Y_tfm)
  sigma = sqrt(diag(Y_cov))
  # regularize correlation matrix
  Y_cor = cor(Y_tfm)
  # cholesky decomposition of correlation matrix
  L = t(chol(Y_cor))

  # covariance matrix across donors
  Y_donor = df_samples_subset %>%
    group_by_at(group) %>%
    summarise_at(protein_names, median) %>%
    dplyr::select(protein_names)
  Y_donor_svd = Y_donor %>% tfm %>% cov %>% svd
  #r_donor = sum(cumsum(Y_donor_svd$d/sum(Y_donor_svd$d)) < 0.95)
  sigma_donor = sqrt(Y_donor_svd$d[1:r_donor])
  Q_donor = Y_donor_svd$u[,1:r_donor]
  x_donor = c(Q_donor)

  # sample zeroinflation estimate
  theta = df_samples_subset %>%
    group_by_(group) %>%
    summarize_at(protein_names, function(x) mean(x == 0)) %>%
    dplyr::select(protein_names) %>%
    as.matrix %>%
    t
  theta[which(theta == 0, arr.ind = TRUE)] = 0.001

  # set random effects to zero
  z = rep(0, n*d);
  z_donor = rep(0, k*r_donor)
  stan_init = list(
    beta = beta,
    # cell level
    sigma = sigma,
    L = L,
    z = z,
    # donor level
    sigma_donor = sigma_donor,
    x_donor = x_donor,
    z_donor = z_donor,
    # zero inflation
    theta = theta
  )

  # compile model
  stan_file = system.file("exec", "poisson.stan", package = "cytoeffect")
  #stan_file = "../../exec/poisson.stan"
  model = stan_model(file = stan_file, model_name = "poisson")

  # # run sampler
  # fit_mle = optimizing(model,
  #                      data = stan_data,
  #                      init = stan_init,
  #                      as_vector = FALSE,
  #                      verbose = TRUE)
  # stan_init = fit_mle$par[c("beta",
  #                           "sigma","sigma_term","sigma_donor",
  #                           "z","z_term","z_donor",
  #                           "x","x_term","x_donor")]
  fit_mcmc = sampling(model,
                      pars = c("beta",
                               "sigma",
                               "sigma_donor",
                               "Q_donor",
                               "Cor",
                               "Cor_donor",
                               "b_donor",
                               "theta"
                               ),
                      data = stan_data,
                      iter = iter,
                      warmup = warmup,
                      chains = num_chains,
                      cores = num_chains,
                      seed = seed,
                      init = rep(list(stan_init), num_chains),
                      save_warmup = FALSE,
                      control = list(adapt_delta = adapt_delta))

  # # Laplace approximation
  # stan_file = system.file("exec", "poisson_eb.stan", package = "cytoeffect")
  # #stan_file = "../../exec/poisson_eb.stan"
  # model_eb = stan_model(file = stan_file, model_name = "poisson_eb")
  # stan_data = list(Y = Y, n = n, d = d, p = p,
  #                  k = k, donor = donor, term = term,
  #                  r = rank,
  #                  b = fit_mle$par$b)
  # fit_mle = optimizing(model_eb,
  #                      data = stan_data,
  #                      init = stan_init,
  #                      as_vector = FALSE,
  #                      hessian = TRUE,
  #                      verbose = TRUE)
  # neg_hessian = -fit_mle$hessian
  # cond1 = paste0("beta.",1:10,".1")
  # cond2 = paste0("beta.",1:10,".2")
  # beta_names = c(cond1, cond2)
  # beta_sd = sqrt(1/diag(neg_hessian[beta_names, beta_names]))
  # tibble(
  #   protein_names,
  #   beta_sd[cond1],
  #   beta_sd[cond2]
  # )

  # create cytoeffect class
  obj = list(fit_mcmc = fit_mcmc,
             #fit_mle = fit_mle,
             df_samples_subset = df_samples_subset,
             protein_names = protein_names,
             condition = condition,
             group = group)
  class(obj) = "cytoeffect_poisson"
  obj

}
