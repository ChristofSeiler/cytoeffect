context("data simulation")

test_that("stop when simulated data is not stored in a data frame", {

  testthat::expect_is(simulate_data(), "data.frame")

})
