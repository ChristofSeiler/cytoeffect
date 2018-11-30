#' HMC sampling for Poisson Log-Normal Model.
#'
#' @import rstan
#' @import reshape2
#' @import dplyr
#' @import batchtools
#' @export
#'
poisson_lognormal = function(df_samples_subset,
                             protein_names,
                             condition,
                             group,
                             iter = 325,
                             warmup = 200,
                             num_chains = 4,
                             eta = 1.0) {

  # some checks
  if(sum(names(df_samples_subset) == condition) == 0)
    stop("condition column missing")
  if(sum(names(df_samples_subset) == group) == 0)
    stop("group column missing")
  if(nrow(df_samples_subset) == 0)
    stop("no observations")
  if(!eta > 0)
    stop("eta needs to be positive")

  # prepare input data
  #df_samples_subset$group_condition = paste0(pull(df_samples_subset, group), "_",
  #                                      pull(df_samples_subset, condition))
  Y = df_samples_subset %>% dplyr::select(protein_names) %>% as.matrix()
  #X = model.matrix(formula(paste("~",condition,"*celltype")), data = df_samples_subset)
  X = model.matrix(formula(paste("~",condition)), data = df_samples_subset)
  n = nrow(Y)
  d = ncol(Y)
  p = ncol(X)
  #k = length(unique(df_samples_subset$group_condition))
  #donor = as.integer(as.factor(df_samples_subset$group_condition))
  donor = df_samples_subset %>%
    pull(group) %>%
    as.factor() %>%
    as.integer()
  k = length(table(donor))
  stan_data = list(Y = Y, X = X, n = n, d = d, p = p,
                   k = k, donor = donor, eta = eta)

  # prepare starting point for sampler
  beta = lapply(1:ncol(Y), function(j) {
    coefs = glm.fit(X, Y[,j], family = stats::poisson())$coefficients
    tibble(x0 = coefs[1], x1 = coefs[2])
  }) %>% bind_rows()
  logY = log(Y + 1)
  sigma = sqrt(diag(cov(logY)))
  L = t(chol(cor(logY)))
  z = matrix(0, nrow = n, ncol = d)
  z_donor = matrix(0, nrow = k, ncol = d)
  stan_init = list(beta = beta,
                   sigma = sigma, sigma_donor = sigma,
                   L = L, L_donor = L,
                   z = z, z_donor = z_donor)

  # cluster function
  run_sampling = function(seed) {

    # compile model
    stan_file = system.file("exec", "poisson.stan", package = "cytoeffect")
    model = stan_model(file = stan_file, model_name = "poisson")

    # run sampler
    fit_mcmc = sampling(model,
                        pars = c("beta",
                                 "sigma","sigma_term","sigma_donor",
                                 "L","L_term","L_donor",
                                 "Cor","Cor_term","Cor_donor",
                                 "b_donor",
                                 "Y_hat"
                                 ),
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = num_chains,
                        cores = num_chains,
                        seed = seed,
                        init = rep(list(stan_init), num_chains))
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
             protein_names = protein_names,
             conditions = levels(pull(df_samples_subset, condition)),
             celltypes = levels(as.factor(df_samples_subset$celltype)),
             covariates = colnames(X),
             Y = Y)
  class(obj) = "cytoeffect_poisson"
  obj

}
