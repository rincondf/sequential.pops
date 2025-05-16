test_that("Evaluation of STBP works", {
  counts3 <- rpois(20, lambda = 3)

  test2F <- stbp_composite(data = counts3,
                           greater_than = TRUE,
                           hypothesis = 5,
                           density_func = "poisson",
                           prior = 0.5,
                           lower_bnd = 0,
                           upper_bnd = Inf,
                           lower_criterion = 0.01,
                           upper_criterion = 0.99)
  range <- seq(1, 8)

  expect_length(STBP.eval(test2F,
                          eval.range = seq(1, 8),
                          n = 1,
                          prior = 0.5,
                          N = 30)$AvgSamples, length(range))
})
