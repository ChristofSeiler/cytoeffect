context("poisson model fit")

test_that("stop when input df_sample_subset has no observations", {

  df_samples_subset = tibble(
    M1 = integer(),
    M2 = integer(),
    condition = factor(),
    group = factor()
    )
  protein_names = c("M1","M2")
  condition = "condition"
  group = "group"
  ncores = 1

  expect_error(
    cytoeffect::poisson_lognormal(df_samples_subset, protein_names,
                                  condition = condition, group = group,
                                  num_chains = ncores)
  )

})

test_that("stop when input df_sample_subset has no condition column", {

  df_samples_subset = tibble(
    M1 = integer(),
    M2 = integer(),
    group = factor()
  )
  protein_names = c("M1","M2")
  condition = "condition"
  group = "group"
  ncores = 1

  expect_error(
    cytoeffect::poisson_lognormal(df_samples_subset, protein_names,
                                  condition = condition, group = group,
                                  num_chains = ncores)
  )

})

test_that("stop when input df_sample_subset has no group column", {

  df_samples_subset = tibble(
    M1 = integer(),
    M2 = integer(),
    condition = factor()
  )
  protein_names = c("M1","M2")
  condition = "condition"
  group = "group"
  ncores = 1

  expect_error(
    cytoeffect::poisson_lognormal(df_samples_subset, protein_names,
                                  condition = condition, group = group,
                                  num_chains = ncores)
  )

})

test_that("stop when input df_sample_subset has a condition with more than two levels", {

  df_samples_subset = tibble(
    M1 = c(3,2,1),
    M2 = c(1,1,1),
    condition = factor(c("unstim","stim","control")),
    group = factor(c("A","B","B"))
  )
  protein_names = c("M1","M2")
  condition = "condition"
  group = "group"
  ncores = 1

  expect_error(
    cytoeffect::poisson_lognormal(df_samples_subset, protein_names,
                                  condition = condition, group = group,
                                  num_chains = ncores)
  )

})

test_that("fit poisson model", {

  set.seed(0xdada)
  latent_M = rpois(100, 3)
  df_samples_subset = tibble(
    M1 = rpois(100, latent_M + 5),
    M2 = rpois(100, latent_M + 4),
    condition = factor(c(rep("unstim",50), rep("stim",50))),
    group = factor(c(rep("A",25), rep("B",25), rep("C",25), rep("D",25)))
  )
  protein_names = c("M1","M2")
  condition = "condition"
  group = "group"
  ncores = 1

  obj = cytoeffect::poisson_lognormal(df_samples_subset, protein_names,
                                      condition = condition, group = group,
                                      num_chains = ncores)

  expect_is(obj$fit_mcmc, "stanfit")

})
