% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STBP.R
\name{stbp_simple}
\alias{stbp_simple}
\title{Sequential test of Bayesian posterior probabilities for simple hypotheses
about species absence}
\usage{
stbp_simple(
  data,
  density_func,
  overdispersion = NA,
  prior = 0.5,
  upper_bnd = Inf,
  lower_criterion,
  upper_criterion
)
}
\arguments{
\item{data}{For count data, either a vector (for purely sequential designs) o
a matrix (group sequential designs) with sequential count data, with sampling
bouts collected over time in columns and samples within bouts in rows. NAs
are allowed in case sample size within bouts is unbalanced. For binomial
data, a list of matrices with integer non-negative values of observations
in col 1 and number of samples in col 2, so that each matrix within the list
corresponds to a sampling bout. NAs are \emph{not} allowed for binomial data.}

\item{density_func}{Kernel probability density function for the data. See details.}

\item{overdispersion}{A character string (if a function) or a number
specifying the overdispersion parameter. Only required when using
\code{"negative binomial"} or \code{"beta-binomial"} as kernel densities.
See details.}

\item{prior}{Single number with initial prior. Must be in the interval
\eqn{[0,1]}. If no prior information is available 0.5 (default) is recommended.}

\item{upper_bnd}{Single number indicating the greatest possible value for \eqn{\mu}.
For count data, is often \code{Inf} (default), but it must be \eqn{\leq 1}
for binomial data.}

\item{lower_criterion}{Criterion to decide against the tested hypothesis.
This is the lowest credibility to the hypothesis to stop sampling and decide
against.}

\item{upper_criterion}{Criterion to decide in favor of the tested hypothesis.
This is the greatest credibility to the hypothesis to stop sampling and decide
in favor.}
}
\value{
An object of class \code{"STBP"}.
}
\description{
Runs a Sequential test of Bayesian Posterior Probabilities for hypotheses
about species absence of the form \eqn{H:\mu = 0}. Data is treated in a
sequential framework.
}
\details{
The \code{density_func} argument should be specified as character string.
Acceptable options are \code{"poisson"}, \code{"negative binomial"},
\code{"binomial"} and \code{"beta-binomial"}. The overdispersion
parameter for \code{"negative binomial"} and \code{"beta-binomial"} can be
either a constant or a function of the mean.

If a function, it should be specified as a character string with the name of
an existing function. For options of empirical functions to describe
overdispersion as a function of the mean see Binns et al. (2000). The most
common approach for the negative binomial family is Taylor's Power Law, which
describes the variance as a function of the mean with two parameters, \eqn{a}
and \eqn{b}. Overdispersion, \eqn{k}, can then be specified as: \deqn{k = \frac{\mu^2}{a \mu^b - \mu}}
}
\examples{
# Testing the absence of a species in a given area from a sequential random
# sampling of 3 bouts made of 10 samples (counts) each (all absences). Upper
# criterion set to 0.9999

counts10 <- matrix(rep(0, 30), 10, 3)

test1G <- stbp_simple(data = counts10,
                        density_func = "poisson",
                        prior = 0.5,
                        upper_bnd = Inf,
                        lower_criterion = 0,
                        upper_criterion = 0.9999)
test1G

# returns a recommendation of "keep sampling" due to insufficient evidence.

# Testing the same hypothesis with the same upper criterion but from a
# sequential random sampling of 3 bouts made of 30 samples (counts) each
# (all absences).

counts30 <- matrix(rep(0, 90), 30, 3)

test2G <- stbp_simple(data = counts30,
                        density_func= "poisson",
                        prior = 0.5,
                        upper_bnd = Inf,
                        lower_criterion = 0,
                        upper_criterion = 0.9999)
test2G

# returns a recommendation of "accept H" of the species being absent from
# that area.

## End (Not run)
}
\references{
Binns, M.R., Nyrop, J.P. & Werf, W.v.d. (2000) \emph{Sampling and
monitoring in crop protection: the theoretical basis for developing practical
decision guides}. CABI Pub., Wallingford, Oxon, UK; New York, N.Y.

Rincon, D.F., McCabe, I. & Crowder, D.W. (2025) Sequential testing of
complementary hypotheses about population density. \emph{Methods in Ecology
and Evolution}. <https://doi.org/10.1111/2041-210X.70053>
}
