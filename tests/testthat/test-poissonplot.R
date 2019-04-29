context("poisson model plots")

test_that("plot beta from poisson model", {

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

  expect_is(plot(obj, type = "beta"), "ggplot")
  expect_is(plot(obj, type = "sigma"), "ggplot")
  expect_is(plot(obj, type = "Cor")[[1]], "ggplot")
  expect_is(plot(obj, type = "Cor")[[2]], "ggplot")
  expect_is(plot(obj, type = "Cor")[[3]], "ggplot")
  expect_error(plot(obj, type = "something_crazy"))

  expect_is(plot_mds(obj, asp = TRUE, ncores = 1), "ggplot")
  expect_is(plot_mds(obj, asp = FALSE, ncores = 1), "ggplot")
  expect_error(plot_mds(obj$fit_mcmc))

})
