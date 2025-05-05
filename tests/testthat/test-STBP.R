test_that("Posterior calculation for composites works", {
  counts <- c(1, 2, 3)

  expect_equal(stbp_posterior_composite(data = counts,
                                        greater_than = TRUE,
                                        hypothesis = 2,
                                        likelihood_func = function(data, x)
                                            {dpois(data, lambda = x)},
                                        prior = 0.5,
                                        lower_bnd = 0,
                                        upper_bnd = Inf), 0.60630278)
})


test_that("Posterior calculation for composites works", {
  counts <- c(1, 2, 3)

  expect_equal(stbp_posterior_composite(data = counts,
                                        greater_than = TRUE,
                                        hypothesis = 2,
                                        likelihood_func = function(data, x)
                                            {dnbinom(data, size = 2, mu = x)},
                                        prior = 0.5,
                                        lower_bnd = 0,
                                        upper_bnd = Inf), 0.72558593)
})


test_that("STBP for composites works", {
  set.seed(101)
  counts3 <- rpois(5, lambda = 3)

  expect_equal(stbp_composite(data = counts3,
                              greater_than = TRUE,
                              hypothesis = 5,
                              likelihood_func = function(data, x)
                                  {dpois(data, lambda = x)},
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = Inf,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)$recommendation,
               "reject H")
})



test_that("STBP for composites works", {
  set.seed(101)
  counts3 <- rpois(5, lambda = 3)

  expect_equal(stbp_composite(data = counts3,
                              greater_than = TRUE,
                              hypothesis = 5,
                              likelihood_func = function(data, x)
                                  {dpois(data, lambda = x)},
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = Inf,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)$num_iterations,
               2)
})





test_that("STBP for composites works", {

  H <- c(2, 5, 10, 20, 40, 40, 20, 10, 5, 2)
  countP <- matrix(NA, 3, 10)

  set.seed(101)
  for(i in 1:10){
    countP[, i] <- rpois(3, lambda = (H[i] - 1))
  }

  expect_equal(stbp_composite(data = countP,
                              greater_than = TRUE,
                              hypothesis = H,
                              likelihood_func = function(data, x)
                                  {dpois(data, lambda = x)},
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = Inf,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)$recommendation,
               "reject H")
})


test_that("STBP for composites works", {

  H <- c(2, 5, 10, 20, 40, 40, 20, 10, 5, 2)
  countP <- matrix(NA, 3, 10)

  set.seed(101)
  for(i in 1:10){
    countP[, i] <- rpois(3, lambda = (H[i] - 1))
  }

  expect_equal(stbp_composite(data = countP,
                              greater_than = TRUE,
                              hypothesis = H,
                              likelihood_func = function(data, x)
                                  {dpois(data, lambda = x)},
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = Inf,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)$num_iterations,
               8)
})


test_that("Posterior calculation for simples works", {
  counts <- c(0, 0, 0)

  expect_equal(stbp_posterior_simple(data = counts,
                                     hypothesis = 0,
                                     likelihood_func = function(data, x)
                                     {dpois(data, lambda = x)},
                                     prior = 0.5,
                                     upper_bnd = Inf), 0.75)
})




test_that("STBP for simples works", {
  counts10 <- matrix(rep(0, 30), 10, 3)

  expect_equal(stbp_simple(data = counts10,
                           hypothesis = 0,
                           likelihood_func= function(data, x)
                           {dpois(data, lambda = x)},
                           prior = 0.5,
                           upper_bnd = Inf,
                           lower_criterion = 0,
                           upper_criterion = 0.9999)$recommendation,
               "keep sampling")
})


test_that("STBP for simples works", {
  counts10 <- matrix(rep(0, 30), 10, 3)

  expect_equal(stbp_simple(data = counts10,
                           hypothesis = 0,
                           likelihood_func= function(data, x)
                           {dpois(data, lambda = x)},
                           prior = 0.5,
                           upper_bnd = Inf,
                           lower_criterion = 0,
                           upper_criterion = 0.9999)$num_iterations,
               3)
})







test_that("STBP for simples works", {

  counts30 <- matrix(rep(0, 90), 30, 3)

  expect_equal(stbp_simple(data = counts30,
                           hypothesis = 0,
                           likelihood_func= function(data, x)
                           {dpois(data, lambda = x)},
                           prior = 0.5,
                           upper_bnd = Inf,
                           lower_criterion = 0,
                           upper_criterion = 0.9999)$recommendation,
               "accept H")
})


test_that("STBP for simples works", {

  counts30 <- matrix(rep(0, 90), 30, 3)

  expect_equal(stbp_simple(data = counts30,
                           hypothesis = 0,
                           likelihood_func= function(data, x)
                           {dpois(data, lambda = x)},
                           prior = 0.5,
                           upper_bnd = Inf,
                           lower_criterion = 0,
                           upper_criterion = 0.9999)$num_iterations,
               3)
})








