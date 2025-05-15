test_that("Sequential probability ratio test works", {

  counts <- c(2, 5, 6, 2, 7)

  expect_equal(sprt(data = counts,
                    mu0 = 2,
                    mu1 = 4,
                    density_func = "negative binomial",
                    overdispersion = 4.6,
                    alpha = 0.1,
                    beta = 0.1)@recommendation, "accept H1")
})

test_that("Sequential probability ratio test works2", {

  counts <- c(2, 5, 6, 2, 7)

  expect_equal(sprt(data = counts,
                    mu0 = 2,
                    mu1 = 4,
                    density_func = "negative binomial",
                    overdispersion = 4.6,
                    alpha = 0.1,
                    beta = 0.1)@iterations, 5)
})
