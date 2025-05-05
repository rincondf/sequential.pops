#' Posterior probability calculation for composite hypotheses
#'
#' @description
#' This function calculates a posterior probability for hypotheses about population densities of the form \eqn{H:\mu > \psi} or \eqn{H:\mu < \psi}, given the data at a single
#' iteration. This function is to be used in a sequential framework, and called on the sequential test \code{\link{stbp_composite}}.
#'
#'
#' @param data Vector with count data for a single sampling bout.
#' @param greater_than logical; if TRUE (default), the tested hypothesis is of the form \eqn{H:\mu > \psi} otherwise, \eqn{H:\mu < \psi}.
#' @param hypothesis Single value with the hypothesized value of \eqn{\mu}.
#' @param likelihood_func Kernel probability density function for the data. See details.
#' @param prior Single number with initial prior. Must be in the interval \eqn{[0,1]}
#' @param lower_bnd Single number indicating the lower bound of the parameter space for \eqn{\mu}. Most cases is \eqn{0} (default).
#' @param upper_bnd Single number indicating the upper bound of the parameter space for \eqn{\mu}. Most cases is \code{Inf} (default).
#' @returns A single probability value
#' @details
#' The \code{likelihood_func} argument should be specified as a density function with an \dQuote{x} in the
#'  mean argument and \dQuote{data} in the vector of quantiles.
#'  For example, \code{function(data, x){dpois(data, lambda = x)}}.
#' @examples
#'
#' # Counts collected in a single sampling bout
#' counts <- c(1, 2, 3)
#'
#' # Calculate posterior probability from a naive 0.5 prior for H1:mu>2
#' # (a population being >2 individuals per sampling unit) with
#' # a poisson kernel
#'
#' stbp_posterior_composite(data = counts,
#'                           greater_than = TRUE,
#'                           hypothesis = 2,
#'                           likelihood_func = function(data, x)
#'                                                 {dpois(data, lambda = x)},
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf) # returns 0.60630278
#'
#' # Same analysis but with a negative binomial kernel.
#' # Note that 'size' can either be a positive number or a function.
#'
#' stbp_posterior_composite(data = counts,
#'                           greater_than = TRUE,
#'                           hypothesis = 2,
#'                           likelihood_func = function(data, x)
#'                                                 {dnbinom(data, size = 2,
#'                                                          mu = x)},
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf) # returns 0.72558593
#' @export

stbp_posterior_composite <- function(data,
                                     greater_than = TRUE,
                                     hypothesis,
                                     likelihood_func,
                                     prior,
                                     lower_bnd = 0,
                                     upper_bnd = Inf) {
  likelihood <- function(x) {
    prod(likelihood_func(data, x))
  }

  H1 <- prior *
    integrate(
      Vectorize(likelihood),
      lower = if(greater_than == TRUE) hypothesis else lower_bnd,
      upper = if(greater_than == TRUE) upper_bnd else hypothesis
    )$value
  H0 <- (1 - prior) *
    integrate(
      Vectorize(likelihood),
      lower = if(greater_than == TRUE) lower_bnd else hypothesis,
      upper = if(greater_than == TRUE) hypothesis else upper_bnd
    )$value
  posterior <- H1 / (H0 + H1)
  posterior
}

#' Sequential test of Bayesian posterior probabilities for composite hypotheses
#'
#' @description
#' Runs a Sequential test of Bayesian Posterior Probabilities for hypotheses
#' about population densities of the form \eqn{H:\mu > \psi} or \eqn{H:\mu < \psi}.
#' Data is treated in a sequential framework.
#'
#'
#' @param data Either a vector (for purely sequential designs) o a matrix
#' (group sequential designs) with sequential count data, with sampling bouts
#' collected over time in columns and sampling within bouts in rows.
#' @param greater_than logical; if TRUE (default), the tested hypothesis is of the form \eqn{H:\mu > \psi} otherwise, \eqn{H:\mu < \psi}.
#' @param hypothesis Either a single value or a vector with the hypothesized values for \eqn{\mu}.
#' If a vector, should contain at least as many values as \code{ncol(data)}
#' @param likelihood_func Kernel probability density function for the data. See details.
#' @param prior Single number with initial prior. Must be in the interval \eqn{[0,1]}.
#' @param lower_bnd Single number indicating the lower bound of the parameter space for \eqn{\mu}. Most cases is \eqn{0} (default).
#' @param upper_bnd Single number indicating the upper bound of the parameter space for \eqn{\mu}. Most cases is \code{Inf} (default).
#' @param lower_criterion Criterion to decide against the tested hypothesis.
#' This is the lowest credibility to the hypothesis to stop sampling and decide against.
#' @param upper_criterion Criterion to decide in favor of the tested hypothesis.
#' This is the greatest credibility to the hypothesis to stop sampling and decide in favor.
#'
#' @details
#' The \code{likelihood_func} argument should be specified as a density function
#' with an \dQuote{x} in the mean argument and \dQuote{data} in the vector of quantiles.
#' For example, \code{function(data, x){dpois(data, lambda = x)}}.
#'
#' @returns
#' A list with a vector of posterior \code{$probabilities} of length
#' \code{ncol(data)}, \code{num_iterations} with the number of iterations
#' required to reach a decision, and a \code{$recommendation}, either
#' \dQuote{reject H}, \dQuote{accept H}, or \dQuote{keep sampling} for a decision in favor,
#' against or inconclusive about the hypothesis. The test stops whenever the posterior
#' probability reached either of the upper or lower criteria.
#'
#' @examples
#' # Testing the hypothesis of a population size being greater than 5 individuals
#' # per sampling unit (H: mu > 5). The sequential sampling is made of 5 sampling
#' # bouts made of one sample each.
#'
#' set.seed(101)
#' counts3 <- rpois(5, lambda = 3)
#'
#' stbp_composite(data = counts3,
#'                 greater_than = TRUE,
#'                 hypothesis = 5,
#'                 likelihood_func = function(data, x)
#'                     {dpois(data, lambda = x)},
#'                 prior = 0.5,
#'                 lower_bnd = 0,
#'                 upper_bnd = Inf,
#'                 lower_criterion = 0.001,
#'                 upper_criterion = 0.999)$recommendation
#' # returns "reject H" after 2 iterations.
#'
#' # Testing the hypothesis of a sampled population being greater than trajectory H
#' H <- c(2, 5, 10, 20, 40, 40, 20, 10, 5, 2)
#'
#' # Generating sequential samples (n = 3) from a population that is 1 below H
#' # (H - 1)
#'
#' countP <- matrix(NA, 3, 10)
#' set.seed(101)
#' for(i in 1:10){
#'   countP[, i] <- rpois(3, lambda = (H[i] - 1))
#' }
#'
#' # Running STBP on the sample
#'
#' stbp_composite(data = countP,
#'                 greater_than = TRUE,
#'                 hypothesis = H,
#'                 likelihood_func = function(data, x)
#'                     {dpois(data, lambda = x)},
#'                 prior = 0.5,
#'                 lower_bnd = 0,
#'                 upper_bnd = Inf,
#'                 lower_criterion = 0.001,
#'                 upper_criterion = 0.999)$recommendation
#' # returns "reject H" after 8 iterations.
#' @export
#'
stbp_composite <- function(data,
                           greater_than = TRUE,
                           hypothesis,
                           likelihood_func,
                           prior = 0.5,
                           lower_bnd = 0,
                           upper_bnd = Inf,
                           lower_criterion = 0.01,
                           upper_criterion = 0.99) {

  # useful to treat data as a matrix to to able to process group or
  # single sequential data

  if(is.vector(data)) data <- matrix(data, 1, length(data))

  # If hypothesis is just a single repeated value,
  # make a vector of that value repeated as many times as there are bouts.
  # This makes it so that the user can input either a hypothesis made of
  # a single value or a trajectory.

  if(length(hypothesis) == 1) hypothesis <- rep(hypothesis, ncol(data))

  # Init vector with length equal to number of sampling bouts
  # and with initial prior as its first value
  posteriors <- c(prior, rep(NA, ncol(data) - 1))
  for(i in 1: ncol(data)) {
    if(posteriors[i] < 0.001) posteriors[i] <- 0.001
    if(posteriors[i] > 0.999) posteriors[i] <- 0.999
    bout = data[, i]
    posteriors[i + 1] = stbp_posterior_composite(bout,
                                                 greater_than,
                                                 hypothesis[i],
                                                 likelihood_func,
                                                 prior = posteriors[i],
                                                 lower_bnd = lower_bnd,
                                                 upper_bnd = upper_bnd)
    # Break from iteration early if early_return is true,
    # the minimum iterations have been reached, and
    # if either of the decision criteria have been reached
    if(((posteriors[i + 1] < lower_criterion) ||
        (posteriors[i + 1] > upper_criterion))
    ) break
  }

  if(((posteriors[i + 1] < lower_criterion) ||
      (posteriors[i + 1] > upper_criterion))
  ) {
    if(posteriors[i + 1] < lower_criterion) response <- "reject H"
    else response <- "accept H"
  }
  else response <- "keep sampling"

  return(list(
    probabilities = posteriors,
    recommendation = response,
    num_iterations = i)
  )
}


#' Posterior calculation for simple hypotheses
#'
#' @description
#' This function calculates a posterior probability for hypotheses about population
#' densities, of the form \eqn{H:\mu = \psi}, given the data at a single
#' iteration. This function is to be used in a sequential framework, and called
#' on the sequential test \code{\link{stbp_simple}}.
#'
#' @param data Vector with count data for a single sampling bout.
#' @param hypothesis Single value with the hypothesized value of \eqn{\mu}.
#' @param likelihood_func Kernel probability density function for the data. See details.
#' @param prior Single number with initial prior. Must be in the interval \eqn{[0,1]}
#' @param upper_bnd Single number indicating the greatest possible value for \eqn{\mu}. Most cases is \code{Inf}.
#'
#' @return A single probability value
#'
#' @details
#' The \code{likelihood_func} argument should be specified as a density function with an \dQuote{x} in the
#'  mean argument and \dQuote{data} in the vector of quantiles.
#'  For example, \code{function(data, x){dpois(data, lambda = x)}}.
#' @examples
#' # Counts collected in a single sampling bout
#' counts <- c(0, 0, 0)
#'
#' # Calculate posterior probability from a naive 0.5 prior for H:mu=0
#' # (a species being absent in an area) with a poisson kernel.
#'
#' stbp_posterior_simple(data = counts,
#'                        hypothesis = 0,
#'                        likelihood_func = function(data, x)
#'                                                 {dpois(data, lambda = x)},
#'                        prior = 0.5,
#'                        upper_bnd = Inf) # returns 0.75
#' @export

stbp_posterior_simple <- function(data,
                                 hypothesis,
                                 likelihood_func,
                                 prior,
                                 upper_bnd = Inf) {
  likelihood <- function(x) {
    prod(likelihood_func(data, x))
  }

  H1 <- prior *
    likelihood(hypothesis)
  H0 <- (1 - prior) *
    integrate(
      Vectorize(likelihood),
      lower = hypothesis,
      upper = Inf
    )$value
  posterior <- H1 / (H0 + H1)
  posterior
}

#' Sequential test of Bayesian posterior probabilities for simple hypotheses
#'
#' @description
#' Runs a Sequential test of Bayesian Posterior Probabilities for hypotheses about population densities of the form \eqn{H:\mu = \psi}.
#' Data is treated in a sequential framework.
#'
#'
#' @param data Either a vector (for purely sequential designs) o a matrix
#' (group sequential designs) with sequential count data, with sampling bouts
#' collected over time in columns and sampling within bouts in rows.
#' @param hypothesis Either a single value or a vector with the hypothesized values for \eqn{\mu}.
#' If a vector, should contain at least as many values as \code{ncol(data)}
#' @param likelihood_func Kernel probability density function for the data. See details.
#' @param prior Single number with initial prior. Must be in the interval \eqn{[0,1]}.
#' @param upper_bnd Single number indicating the greatest possible value for \eqn{\mu}. Most cases is \code{Inf} (default).
#' @param lower_criterion Criterion to decide against the tested hypothesis.
#' This is the lowest credibility to the hypothesis to stop sampling and decide against.
#' @param upper_criterion Criterion to decide in favor of the tested hypothesis.
#' This is the greatest credibility to the hypothesis to stop sampling and decide in favor.
#'
#' @details
#' The \code{likelihood_func} argument should be specified as a density function
#' with an \dQuote{x} in the mean argument and \dQuote{data} in the vector of quantiles.
#' For example, \code{function(data, x){dpois(data, lambda = x)}}.
#'
#' @returns
#'
#' A list with a vector of posterior \code{$probabilities} of length
#' \code{ncol(data)}, \code{num_iterations} with the number of iterations
#' required to reach a decision, and a \code{$recommendation}, either
#' \dQuote{reject H}, \dQuote{accept H}, or \dQuote{keep sampling} for a decision in favor,
#' against or inconclusive about the hypothesis. The test stops whenever the posterior
#' probability reached either of the upper or lower criteria.
#'
#' @examples
#'
#' # Testing the absence of a species in a given area from a sequential random
#' # sampling of 3 bouts made of 10 samples (counts) each (all absences). Upper
#' # criterion set to 0.9999
#'
#' counts10 <- matrix(rep(0, 30), 10, 3)
#'
#' stbp_simple(data = counts10,
#'               hypothesis = 0,
#'               likelihood_func= function(data, x)
#'                 {dpois(data, lambda = x)},
#'               prior = 0.5,
#'               upper_bnd = Inf,
#'               lower_criterion = 0,
#'               upper_criterion = 0.9999)$recommendation
#'
#' # returns a recommendation of "keep sampling" due to insufficient evidence
#' # after processing the 3 sampling bouts.
#'
#' # Testing the same hypothesis with the same upper criterion but from a
#' # sequential random sampling of 3 bouts made of 30 samples (counts) each
#' # (all absences).
#'
#' counts30 <- matrix(rep(0, 90), 30, 3)
#'
#' stbp_simple(data = counts30,
#'               hypothesis = 0,
#'               likelihood_func= function(data, x)
#'                 {dpois(data, lambda = x)},
#'               prior = 0.5,
#'               upper_bnd = Inf,
#'               lower_criterion = 0,
#'               upper_criterion = 0.9999)$recommendation
#'
#' # returns a recommendation of "accept H" of the species being absent from
#' # that area after processing the 3 sampling bouts.
#'
#'
#' @export

stbp_simple <- function(data,
                         hypothesis,
                         likelihood_func,
                         prior = 0.5,
                         upper_bnd = Inf,
                         lower_criterion,
                         upper_criterion) {

  # useful to treat data as a matrix to to able to process group or
  # single sequential data

  if(is.vector(data)) data <- matrix(data, 1, length(data))

  # If hypothesis is just a single repeated value,
  # make a vector of that value repeated as many times as there are bouts.
  # This makes it so that the user can input either a hypothesis made of a
  # single value or a trajectory.

  if(length(hypothesis) == 1) hypothesis <- rep(hypothesis, ncol(data))

  # Init vector with length equal to number of sampling bouts
  # and with initial prior as its first value
  posteriors <- c(prior, rep(NA, ncol(data) - 1))
  for(i in 1: ncol(data)) {
    bout = data[, i]
    posteriors[i + 1] = stbp_posterior_simple(bout,
                                              hypothesis[i],
                                              likelihood_func,
                                              prior = posteriors[i],
                                              upper_bnd = upper_bnd)
    # Break from iteration early if early_return is true,
    # the minimum iterations have been reached, and
    # if either of the decision criteria have been reached
    if(((posteriors[i + 1] < lower_criterion) ||
        (posteriors[i + 1] > upper_criterion))
    ) break
  }

  if(((posteriors[i + 1] < lower_criterion) ||
      (posteriors[i + 1] > upper_criterion))
  ) {
    if(posteriors[i + 1] < lower_criterion) response <- "reject H"
    else response <- "accept H"
  }
  else response <- "keep sampling"

  return(list(
    probabilities = posteriors,
    recommendation = response,
    num_iterations = i)
  )
}
