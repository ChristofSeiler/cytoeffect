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

  # prepare starting point for sampler

  # log mean per condition (same as Poisson regression)
  beta = df_samples_subset %>%
    group_by_(condition) %>%
    summarise_at(protein_names, function(x) log(mean(x))) %>%
    dplyr::select(protein_names) %>%
    t
  # set according to prior
  is.neginf = function(x) x == -Inf
  beta[is.neginf(beta)] = -7 # beta[j] ~ normal(0, 7);

  # covariance matrix per condition
  initcov = function(tfmY) {
    sigma = sqrt(diag(cov(tfmY)))
    corY = cor(tfmY)
    corY[is.na(corY)] = 0
    if(attr(chol(corY, pivot = TRUE), "rank") == nrow(corY)) {
      L = t(chol(corY))
    } else {
      L = diag(1, nrow(corY))
    }
    list(sigma=sigma,L=L)
  }
  tfm = function(x) asinh(x/5)
  cov1 = initcov(tfm(Y[term == 1,]))
  cov2 = initcov(tfm(Y[term == 2,]))
  # set to low value when zero
  cov1$sigma = ifelse(cov1$sigma == 0, 0.01, cov1$sigma)
  cov2$sigma = ifelse(cov2$sigma == 0, 0.01, cov2$sigma)

  # covariance matrix across donors
  Y_donor = df_samples_subset %>%
    group_by_(group) %>%
    summarise_at(protein_names, median) %>%
    dplyr::select(protein_names)
  # low rank estimation of donor covariance matrix
  out = svd(cov(tfm(Y_donor)))
  rank = 2
  eigval = c(out$d[1:rank], rep(0,length(out$d)-rank))
  reg = out$u %*% diag(eigval) %*% t(out$v)
  #ggcorrplot::ggcorrplot(cor(reg)) +
  #  ggtitle(paste0("rank = ", rank))
  sigma_donor = sqrt(diag(reg))
  cor_donor = cor(reg)

  # set random effects to zero
  z = matrix(0, nrow = n, ncol = d)
  z_term = matrix(0, nrow = n, ncol = d)
  z_donor = matrix(0, nrow = k, ncol = d)
  stan_init = list(
    beta = beta,
    sigma = cov1$sigma, sigma_term = cov2$sigma, sigma_donor = sigma_donor,
    L = cov1$L, L_term = cov2$L,
    z = z, z_term = z_term, z_donor = z_donor
  )
  stan_data = list(Y = Y, n = n, d = d, p = p,
                   k = k, donor = donor, term = term,
                   cor_donor = cor_donor
  )

  # cluster function
  run_sampling = function(seed) {

    # compile model
    stan_file = system.file("exec", "poisson.stan", package = "cytoeffect")
    model = stan_model(file = stan_file, model_name = "poisson")

    # run sampler
    fit_mcmc = sampling(model,
                        pars = c("beta",
                                 "sigma","sigma_term","sigma_donor",
                                 # "L","L_term","L_donor",
                                 "Cor","Cor_term",#"Cor_donor",
                                 "b_donor"
                                 # "Y_hat"
                                 ),
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = num_chains,
                        cores = num_chains,
                        seed = seed,
                        init = rep(list(stan_init), num_chains),
                        save_warmup = FALSE)
    fit_mcmc
  }

  # preapte and submit cluster job
  current_time = Sys.time() %>%
    str_replace_all(":","") %>%
    str_replace_all("-| ","_")
  reg = makeRegistry(file.dir = paste0("registry_",current_time),
                     packages = "rstan")
  batchMap(run_sampling, seed = 1)
  submitJobs()
  waitForJobs()
  fit_mcmc = reduceResultsList()[[1]]

  # create cytoeffect class
  obj = list(fit_mcmc = fit_mcmc,
             df_samples_subset = df_samples_subset,
             protein_names = protein_names,
             condition = condition,
             group = group)
  class(obj) = "cytoeffect_poisson"
  obj

}
