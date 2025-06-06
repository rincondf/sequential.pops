% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPRT-methods.R
\docType{methods}
\name{plot,SPRT,missing-method}
\alias{plot,SPRT,missing-method}
\title{Draws a plot with the cumulative introduced counts contrasted with stop lines
or stop lines when no data is provided}
\usage{
\S4method{plot}{SPRT,missing}(x, y)
}
\arguments{
\item{x}{Created as a result of a call to \code{\link{sprt}}.}

\item{y}{Unused entry.}
}
\value{
When data is provided, a plot with cumulative counts contrasted with stop lines
from a \code{"SPRT"}. When no data is provided, a plot with stop lines.
}
\description{
Method for signature \code{"SPRT"} to display results or stop lines.
}
\examples{

test00 <- sprt(mu0 = 2,
                mu1 = 4,
                density_func = "negative binomial",
                overdispersion = 4.6,
                alpha = 0.1,
                beta = 0.1)
plot(test00) # returns chart with stop lines

counts <- c(2, 5, 6, 2, 7)

test11 <- sprt(data = counts,
                mu0 = 2,
                mu1 = 4,
                density_func = "negative binomial",
                overdispersion = 4.6,
                alpha = 0.1,
                beta = 0.1)
plot(test11) # returns chart with data and stop lines

## End (Not run)
}
