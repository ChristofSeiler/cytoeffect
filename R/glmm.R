#' HMC Sampling for Logistic Linear Mixed Model
#'
#' \code{glmm} uses Hamiltonion Monte Carlo to sample form an extend logistic
#' linear mixed model. We predict the experimental condition from protein
#' marker expressions.
#'
#' @import rstan
#' @import ggplot2
#' @import reshape2
#' @import dplyr
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
#' @param eta Hyperparametr for LKJ prior
#' @return A list of class \emph{cytoeffect} containing
#'   \item{fit_mcmc}{\code{\link[rstan]{rstan}} object}
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
  if(df_samples_subset %>% pull(condition) %>% nlevels != 2)
    stop("condition variables should have two levels")
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
    donor = as.integer(as.factor(pull(df_samples_subset, group))),
    X = X,
    treatment = treatment,
    eta = eta
  )

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
                      seed = 1)

  # create cytoeffect class
  obj = list(fit_mcmc = fit_mcmc,
             protein_names = protein_names,
             conditions = levels(pull(df_samples_subset, condition)),
             X = X)
  class(obj) = "cytoeffect"
  obj

}
