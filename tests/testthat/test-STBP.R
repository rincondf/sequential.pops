test_that("Posterior calculation for composites works", {
  counts <- c(1, 2, 3)

  expect_equal(stbp_posterior_composite(data = counts,
                                        greater_than = TRUE,
                                        hypothesis = 2,
                                        density_func = "poisson",
                                        prior = 0.5,
                                        lower_bnd = 0,
                                        upper_bnd = Inf), 0.60630278)
})


test_that("Posterior calculation for composites works", {
  counts <- c(1, 2, 3)

  expect_equal(stbp_posterior_composite(data = counts,
                                        greater_than = TRUE,
                                        hypothesis = 2,
                                        density_func = "negative binomial",
                                        overdispersion = 2,
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
                              density_func = "poisson",
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = Inf,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)@recommendation,
               "reject H")
})



test_that("STBP for composites works", {
  set.seed(101)
  counts3 <- rpois(5, lambda = 3)

  expect_equal(stbp_composite(data = counts3,
                              greater_than = TRUE,
                              hypothesis = 5,
                              density_func = "poisson",
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = Inf,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)@iterations,
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
                              density_func = "poisson",
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = Inf,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)@recommendation,
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
                              density_func = "poisson",
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = Inf,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)@iterations,
               8)
})




test_that("STBP for composites works (binomial data)", {

  set.seed(101)
  counts4 <- list()

  for(i in 1: 7) {
    counts4[[i]] <- matrix(c(rbinom(5, size = 10, prob = 0.25), rep(10, 5)),
                           5, 2)
  }

  expect_equal(stbp_composite(data = counts4,
                              greater_than = TRUE,
                              hypothesis = 0.2,
                              density_func = "binomial",
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = 1,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)@recommendation,
               "accept H")
})

test_that("STBP for composites works (binomial data)", {

  set.seed(101)
  counts4 <- list()

  for(i in 1: 7) {
    counts4[[i]] <- matrix(c(rbinom(5, size = 10, prob = 0.25), rep(10, 5)),
                           5, 2)
  }

  expect_equal(stbp_composite(data = counts4,
                              greater_than = TRUE,
                              hypothesis = 0.2,
                              density_func = "binomial",
                              prior = 0.5,
                              lower_bnd = 0,
                              upper_bnd = 1,
                              lower_criterion = 0.001,
                              upper_criterion = 0.999)@iterations,
               3)
})





test_that("Posterior calculation for simples works", {
  counts <- c(0, 0, 0)

  expect_equal(stbp_posterior_simple(data = counts,
                                     density_func = "poisson",
                                     prior = 0.5,
                                     upper_bnd = Inf), 0.75)
})




test_that("STBP for simples works", {
  counts10 <- matrix(rep(0, 30), 10, 3)

  expect_equal(stbp_simple(data = counts10,
                           density_func = "poisson",
                           prior = 0.5,
                           upper_bnd = Inf,
                           lower_criterion = 0,
                           upper_criterion = 0.9999)@recommendation,
               "keep sampling")
})


test_that("STBP for simples works", {
  counts10 <- matrix(rep(0, 30), 10, 3)

  expect_equal(stbp_simple(data = counts10,
                           density_func = "poisson",
                           prior = 0.5,
                           upper_bnd = Inf,
                           lower_criterion = 0,
                           upper_criterion = 0.9999)@iterations,
               3)
})







test_that("STBP for simples works", {

  counts30 <- matrix(rep(0, 90), 30, 3)

  expect_equal(stbp_simple(data = counts30,
                           density_func = "poisson",
                           prior = 0.5,
                           upper_bnd = Inf,
                           lower_criterion = 0,
                           upper_criterion = 0.9999)@recommendation,
               "accept H")
})


test_that("STBP for simples works", {

  counts30 <- matrix(rep(0, 90), 30, 3)

  expect_equal(stbp_simple(data = counts30,
                           density_func= "poisson",
                           prior = 0.5,
                           upper_bnd = Inf,
                           lower_criterion = 0,
                           upper_criterion = 0.9999)@iterations,
               3)
})

