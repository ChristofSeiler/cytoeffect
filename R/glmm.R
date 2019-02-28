#' HMC sampling for GLMM with Multivariate Random Effects.
#'
#' @import rstan
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @import batchtools
#' @export
#'
#' @param df_samples_subset Data frame or tibble with proteins counts, cell condition, and group information.
#' @param protein_names A vector of column names of protein to use in the analysis.
#' @param condition The column name of the condition variable.
#' @param group The column name of the group variable.
#' @param iter Number of iteration per chain for the HMC sampler.
#' @param warmup Number of warm up steps per chain for the HMC sampler.
#' @param num_chains Number of HMC chains to run in parallel.
#' @return A list of class \emph{cytoeffect_poisson} containing
#'   \item{fit_mcmc}{Stan fit}
#'   \item{protein_names}{input protein names}
#'   \item{conditions}{input condition variable}
#'   \item{X}{prediction design matrix}
#'
glmm = function(df_samples_subset,
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
                        seed = seed)
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
             X = X)
  class(obj) = "cytoeffect"
  obj

}
