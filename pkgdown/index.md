# sequential.pops: Sequential analysis of biological populations sizes <a href="https://rincondf.github.io/sequential.pops/"><img src="man/figures/logo.png" align="right" height="139" alt="sequential.pops website" /></a>
The `sequential.pops` package provides tools to test sequentially hypotheses about population sizes. Two approaches are included:

1. The sequential test of Bayesian posterior probabilities (Rincon, D.F., McCabe, I. & Crowder, D.W. 2025): This test is based on a sequence of conditional probabilities obtained from the data and can explicitly consider simple or composite hypotheses about static or dynamic population densities and process either purely sequential (one-at-a-time) or group sequential data.
2. The sequential probability ratio test (Wald, A. 1945): This test is based on probability ratios obtained from the likelihood of the data for two contrasted (non-complementary) hypotheses. The SPRT is designed to test composite hypotheses about static populations and process purely sequential data, but there are modifications to make it work for dynamic populations and group sequential data too.

This package includes tools to run and evaluate sequential tests, produce charts with stop lines, and summarize results.

## Useful readings

[A tutorial that illustrates the use `sequential.pops`](https://rincondf.github.io/sequential.pops/articles/Seq-vignette.html)

[Blog post with motivation and some examples](https://methodsblog.com/2025/06/02/simple-intuitive-tool-for-population-management-decision-making/)

[Formulation of the Sequential test of Bayesian posterior probabilities](https://doi.org/10.1111/2041-210X.70053)

[Formulation of the Sequential probability ratio test](https://www.jstor.org/stable/2235829)

<!-- badges: start -->
[![R-CMD-check](https://github.com/rincondf/sequential.pops/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rincondf/sequential.pops/actions/workflows/R-CMD-check.yaml)
[![](https://www.r-pkg.org/badges/version/sequential.pops)](https://cran.r-project.org/package=sequential.pops)
[![](http://cranlogs.r-pkg.org/badges/grand-total/sequential.pops)](https://cran.r-project.org/package=sequential.pops)
<!-- badges: end -->
