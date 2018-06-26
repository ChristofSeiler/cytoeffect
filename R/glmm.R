#' HMC sampling for GLMM with Multivariate Random Effects.
#'
#' @import rstan
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @import batchtools
#' @export
#'
glmm = function(df_samples_subset,
                protein_names,
                condition,
                group,
                iter = 325,
                warmup = 200,
                num_chains = 4,
                eta = 1.0,
                init = NULL) {

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
  protein_names_collapsed = paste(protein_names,collapse = " + ")
  X = unname(model.matrix(formula(paste("~",protein_names_collapsed)),
                          df_samples_subset))
  attr(X, "assign") = NULL
  treatment = as.integer(pull(df_samples_subset, condition))-1
  stan_data = list(
    N = nrow(X),
    P = ncol(X),
    J = length(unique(pull(df_samples_subset, group))),
    K = length(unique(df_samples_subset$celltype)),
    donor = as.integer(as.factor(pull(df_samples_subset, group))),
    celltype = as.integer(as.factor(df_samples_subset$celltype)),
    X = X,
    treatment = treatment,
    eta = eta
  )

  # initialize sampler with previous simulation
  if(class(init) == "cytoeffect") {
    pars = rstan::extract(init$fit_mcmc)
    stan_init = list(
      beta = apply(pars["beta"][[1]],MARGIN = 2,FUN = median),
      Ld = apply(pars["Ld"][[1]],MARGIN = c(2,3),FUN = median),
      sigmad = apply(pars["sigmad"][[1]],MARGIN = 2,FUN = median),
      zd = apply(pars["zd"][[1]],MARGIN = c(2,3),FUN = median),
      Lc = apply(pars["Lc"][[1]],MARGIN = c(2,3),FUN = median),
      sigmac = apply(pars["sigmac"][[1]],MARGIN = 2,FUN = median),
      zc = apply(pars["zc"][[1]],MARGIN = c(2,3),FUN = median),
      ud = apply(pars["ud"][[1]],MARGIN = c(2,3),FUN = median),
      uc = apply(pars["uc"][[1]],MARGIN = c(2,3),FUN = median),
      Cord = apply(pars["Cord"][[1]],MARGIN = c(2,3),FUN = median),
      Corc = apply(pars["Corc"][[1]],MARGIN = c(2,3),FUN = median)
    )
    stan_init = rep(list(stan_init), num_chains)

  } else {
    stan_init = 'random'
  }

  # cluster function
  run_sampling = function(seed) {

    # compile model
    stan_file = system.file("exec", "glmm.stan", package = "cytoeffect")
    model = stan_model(file = stan_file, model_name = "glmm")

    # run sampler
    fit_mcmc = sampling(model,
                        data = stan_data,
                        iter = iter,
                        warmup = warmup,
                        chains = num_chains,
                        cores = num_chains,
                        seed = seed,
                        init = stan_init)
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
             celltypes = levels(as.factor(df_samples_subset$celltype)))
  class(obj) = "cytoeffect"
  obj

}
