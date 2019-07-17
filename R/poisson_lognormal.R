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
#' @param iter Number of iteration per chain for the HMC sampler
#' @param warmup Number of warm up steps per chain for the HMC sampler
#' @param num_chains Number of HMC chains to run in parallel
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
                             r_cell,
                             r_donor,
                             iter = 325,
                             warmup = 200,
                             num_chains = 4) {

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
  #df_samples_subset$group_condition = paste0(pull(df_samples_subset, group), "_",
  #                                      pull(df_samples_subset, condition))
  Y = df_samples_subset %>% dplyr::select(protein_names) %>% as.matrix()
  #X = model.matrix(formula(paste("~",condition)), data = df_samples_subset)
  term = df_samples_subset %>%
    pull(condition) %>%
    as.factor() %>%
    as.integer()
  p = length(table(term))
  n = nrow(Y)
  d = ncol(Y)
  #p = ncol(X)
  #k = length(unique(df_samples_subset$group_condition))
  #donor = as.integer(as.factor(df_samples_subset$group_condition))
  donor = df_samples_subset %>%
    pull(group) %>%
    as.factor() %>%
    as.integer()
  k = length(table(donor))
  stan_data = list(Y = Y, n = n, d = d, p = p,
                   k = k, donor = donor, term = term,
                   r_cell = r_cell, r_donor = r_donor)

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

  # regularize and decompose covariance matrix
  # c is upper bound on condition number:
  # out$d[1]/out$d[length(out$d)]
  # math details: http://lagrange.math.siu.edu/Olive/slch6.pdf
  # c = 100
  # tfm = function(x) asinh(x/5)
  # initcov = function(Y_raw) {
  #   # transform raw counts
  #   Y_tfm = tfm(Y_raw)
  #   # sample standard deviation
  #   Y_cov = cov(Y_tfm)
  #   sigma = sqrt(diag(Y_cov))
  #   # regularize correlation matrix
  #   Y_cor = cor(Y_tfm)
  #   out = svd(Y_cor)
  #   rho = max(0, (out$d[1] - c*out$d[length(out$d)]) / (c - 1) )
  #   Y_cor_reg = 1/(1+rho) * (Y_cor + diag(rho, nrow(Y_cor)))
  #   # cholesky decomposition of correlation matrix
  #   L = t(chol(Y_cor_reg))
  #   list(sigma = sigma, L = L)
  # }

  # covariance matrix across cells per level of condition
  # cov1 = initcov(Y[term == 1,])
  # cov2 = initcov(Y[term == 2,])

  Y_term_svd = svd(cov(Y[term == 1,]))
  sigma = sqrt(Y_term_svd$d[1:r_cell])
  Q = Y_term_svd$u[,1:r_cell]
  x = c(Q)

  Y_term_svd = svd(cov(Y[term == 2,]))
  sigma_term = sqrt(Y_term_svd$d[1:r_cell])
  Q_term = Y_term_svd$u[,1:r_cell]
  x_term = c(Q_term)

  # covariance matrix across donors
  Y_donor = df_samples_subset %>%
    group_by_at(group) %>%
    summarise_at(protein_names, median) %>%
    dplyr::select(protein_names)
  Y_donor_svd = svd(cov(Y_donor))
  sigma_donor = sqrt(Y_donor_svd$d[1:r_donor])
  Q_donor = Y_donor_svd$u[,1:r_donor]
  x_donor = c(Q_donor)

  # ggcorrplot::ggcorrplot(cov1$L %*% t(cov1$L))
  # ggcorrplot::ggcorrplot(cor(tfm(Y[term == 1,])))
  # ggcorrplot::ggcorrplot(cov2$L %*% t(cov2$L))
  # ggcorrplot::ggcorrplot(cor(tfm(Y[term == 2,])))
  # ggcorrplot::ggcorrplot(cor(tfm(Y_donor)))
  # ggcorrplot::ggcorrplot(cov_donor$L %*% t(cov_donor$L))

  # set random effects to zero
  z = matrix(0, nrow = n, ncol = r_cell)
  z_term = matrix(0, nrow = n, ncol = r_cell)
  z_donor = matrix(0, nrow = k, ncol = r_donor)
  stan_init = list(
    beta = beta,
    sigma = sigma, sigma_term = sigma_term, sigma_donor = sigma_donor,
    x = x, x_term = x_term, x_donor = x_donor,
    z = z, z_term = z_term, z_donor = z_donor
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
                               "sigma","sigma_term","sigma_donor",
                               # "L","L_term","L_donor",
                               "Q","Q_term","Q_donor",
                               "Cor","Cor_term","Cor_donor",
                               "b_donor"
                               # "Y_hat"
                               ),
                      data = stan_data,
                      iter = iter,
                      warmup = warmup,
                      chains = num_chains,
                      cores = num_chains,
                      seed = 1,
                      init = rep(list(stan_init), num_chains),
                      save_warmup = FALSE)

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
