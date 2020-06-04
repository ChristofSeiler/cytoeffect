#' Fit simplified Poisson log-normal model
#'
#' \code{poisson_lognormal} uses maximum composite likelihood estimation to fit
#' Poisson log-normal models to each sample.
#'
#' @import dplyr
#' @import batchtools
#' @importFrom magrittr %>% %<>%
#' @importFrom Matrix cov2cor
#' @export
#'
#' @param df_samples_subset Data frame or tibble with proteins counts,
#'   cell condition, and group information
#' @param protein_names A vector of column names of protein to use in the analysis
#' @param condition The column name of the condition variable
#' @param group The column name of the group variable
#' @param ncores Number of CPU cores
#' @param slurm_settings Path to slurm cluster template for \code{\link[batchtools]{batchtools}}
#' @return A list of class \code{cytoeffect_poisson_mcle} containing
#'   \item{tb_args}{output tibble of model fits}
#'   \item{protein_names}{input protein names}
#'   \item{condition}{input condition variable}
#'   \item{group}{input group names}
#'   \item{df_samples_subset}{input df_samples_subset table}
#'
#' @examples
#' set.seed(1)
#' df = simulate_data(n_cells = 10)
#' str(df)
#' fit = poisson_lognormal_mcle(df,
#'                              protein_names = names(df)[3:ncol(df)],
#'                              condition = "condition",
#'                              group = "donor",
#'                              ncores = 1)
#'
poisson_lognormal_mcle = function(df_samples_subset,
                                  protein_names,
                                  condition,
                                  group,
                                  ncores = 1,
                                  slurm_settings = "slurm_batchtools.tmpl") {

  # some checks
  if(sum(names(df_samples_subset) == condition) == 0)
    stop("condition column missing")
  if(sum(names(df_samples_subset) == group) == 0)
    stop("group column missing")
  if(nrow(df_samples_subset) == 0)
    stop("no observations")
  if(df_samples_subset %>% pull(condition) %>% nlevels != 2)
    stop("condition variables should have two levels")

  tb_args = df_samples_subset %>%
    group_by_at(vars(group, condition)) %>%
    tally()

  Y_list = lapply(1:nrow(tb_args), function(i) {
    df_samples_subset %>%
      filter(.data[[group]] == pull(tb_args, group)[i]) %>%
      filter(.data[[condition]] == pull(tb_args, condition)[i]) %>%
      select_at(protein_names) %>%
      as.matrix()
  })
  tb_args %<>% add_column(Y = Y_list)

  current_time = Sys.time() %>%
    str_replace_all(":","") %>%
    str_replace_all("-| ","_")
  reg = makeRegistry(file.dir = paste0("registry_",current_time),
                     conf.file = NA,
                     packages = "cytoeffect")

  if(file.exists(slurm_settings)) {

    # run on cluster using batchtools
    reg$cluster.functions = makeClusterFunctionsSlurm(slurm_settings,
                                                      scheduler.latency = 30,
                                                      fs.latency = 30)
    saveRegistry(reg)
    batchMap(fun = cytoeffect::fit_poilog,
             args = tb_args %>% ungroup %>% dplyr::select(Y),
             more.args = list(ncores = ncores))
    submitJobs(resources = list(ncpus = ncores, # cores per job
                                memory = 4000, # MB memory per job
                                walltime = 360 # minutes per job
    ))
    # expires only after 7 days to handle cluster instability
    waitForJobs(sleep = 60, expire.after = 10080)
    getStatus()
    findErrors()
    getErrorMessages()
    # by setting missing.val to NULL, we impute failed jobs with NULL
    fit_list = reduceResultsList(missing.val = NULL)

    # collect
    tb_args %<>% add_column(fit = fit_list)
    save(tb_args, file = "tb_args.Rdata")

  } else {

    # run locally
    if(.Platform$OS.type == "windows") {
      reg$cluster.functions = makeClusterFunctionsSocket(ncpus = ncores)
    } else {
      reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = ncores)
    }
    batchMap(fun = cytoeffect::fit_poilog,
             args = tb_args %>% ungroup %>% dplyr::select(Y),
             more.args = list(ncores = ncores))
    submitJobs()
    waitForJobs()
    fit_list = reduceResultsList(missing.val = NULL)
    tb_args %<>% add_column(fit = fit_list)

    # cleaning up
    removeRegistry()

  }

  # create cytoeffect_poisson_mcle class
  obj = list(tb_args = tb_args,
             df_samples_subset = df_samples_subset,
             protein_names = protein_names,
             condition = condition,
             group = group)
  class(obj) = "cytoeffect_poisson_mcle"
  obj

}
