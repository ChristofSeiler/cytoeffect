context("poisson model fit")

test_that("stop when input df_sample_subset has no observations", {

  df_samples_subset = tibble(
    A = integer(),
    B = integer(),
    condition = character(),
    group = character()
    )
  protein_names = c("A","B")
  condition = "condition"
  group = "group"
  ncores = 1

  expect_error(
    cytoeffect::poisson_lognormal(df_samples_subset, protein_names,
                                  condition = condition, group = group,
                                  num_chains = ncores)
  )

})
