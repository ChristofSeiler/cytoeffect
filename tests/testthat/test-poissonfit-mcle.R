context("maximum composite likelihood estimation")

test_that("stop when input Y has not enough observations", {

  Y = matrix(0, nrow = 1, ncol = 2)
  expect_error(
    cytoeffect::fit_poilog(Y, ncores = 1)
  )
  Y = matrix(nrow = 0, ncol = 2)
  expect_error(
    cytoeffect::fit_poilog(Y, ncores = 1)
  )

})

test_that("pairwise fit", {

  set.seed(0xdada)
  n_cores = 8
  n = 10  # num of cells
  p = 3   # num of markers
  lambda = rnorm(n*p, 2.3)
  Y = matrix(rpois(n*p, lambda), nrow = n, ncol = p)

  suppressWarnings(
    fit <- cytoeffect::fit_poilog(Y, ncores = 1)
  )

  expect_is(fit, "list")
  expect_true(length(fit$beta) == p)
  expect_true(nrow(fit$Sigma) == p)
  expect_true(ncol(fit$Sigma) == p)

})

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
    cytoeffect::poisson_lognormal_mcle(df_samples_subset, protein_names,
                                       condition = condition, group = group)
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
    cytoeffect::poisson_lognormal_mcle(df_samples_subset, protein_names,
                                       condition = condition, group = group)
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
    cytoeffect::poisson_lognormal_mcle(df_samples_subset, protein_names,
                                       condition = condition, group = group)
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
    cytoeffect::poisson_lognormal_mcle(df_samples_subset, protein_names,
                                       condition = condition, group = group)
  )

})

test_that("fit poisson model", {

  set.seed(0xdada)
  latent_M = rnorm(100, 0, 1)
  df_samples_subset = tibble(
    M1 = rpois(100, exp(latent_M+log(10))),
    M2 = rpois(100, exp(latent_M+log(11))),
    treatment = factor(c(rep("unstim",50), rep("stim",50))),
    patient = factor(c(rep("A",25), rep("B",25), rep("C",25), rep("D",25)))
  )
  protein_names = c("M1","M2")
  condition = "treatment"
  group = "patient"
  ncores = 1

  suppressWarnings(
    fit <- cytoeffect::poisson_lognormal_mcle(df_samples_subset, protein_names,
                                              condition = condition, group = group,
                                              ncores = ncores)
  )

  expect_true(is_tibble(fit$tb_args))

  expect_is(plot(fit, type = "beta"), "ggplot")
  expect_is(plot(fit, type = "sigma"), "ggplot")
  expect_is(plot(fit, type = "Cor"), "ggplot")
  expect_error(plot(fit, type = "something_crazy"))

  expect_is(plot_distatis(fit, ncores = ncores), "ggplot")
  expect_is(plot_distatis(fit, repel = FALSE, ncores = ncores), "ggplot")
  expect_error(plot_distatis(fit$protein_names))

})
