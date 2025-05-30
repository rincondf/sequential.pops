% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STBP-methods.R
\docType{methods}
\name{show,STBP-method}
\alias{show,STBP-method}
\alias{show}
\title{Displays results from a "STBP" object nicely}
\usage{
\S4method{show}{STBP}(object)
}
\arguments{
\item{object}{Created as a result of a call to \code{\link{stbp_simple}} or \code{\link{stbp_composite}}.}
}
\value{
A summary of the test.
}
\description{
Method for signature \code{"STBP"} to show results.
}
\examples{
set.seed(101)
counts3 <- rpois(5, lambda = 3)

test1F <- stbp_composite(data = counts3,
                          greater_than = TRUE,
                          hypothesis = 5,
                          density_func = "poisson",
                          prior = 0.5,
                          lower_bnd = 0,
                          upper_bnd = Inf,
                          lower_criterion = 0.001,
                          upper_criterion = 0.999)
show(test1F)
# returns "reject H".

counts10 <- matrix(rep(0, 30), 10, 3)

test1G <- stbp_simple(data = counts10,
                        density_func= "poisson",
                        prior = 0.5,
                        upper_bnd = Inf,
                        lower_criterion = 0,
                        upper_criterion = 0.9999)
show(test1G)

# returns "keep sampling" due to insufficient evidence.

## End (Not run)
}
