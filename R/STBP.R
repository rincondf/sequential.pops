#' Posterior probability calculation for composite hypotheses
#'
#' @description
#' This function calculates a posterior probability for hypotheses about
#' population densities of the form \eqn{H:\mu > \psi} or \eqn{H:\mu < \psi},
#' given the data at a single iteration. This function is to be used in a
#' sequential framework, and called on the sequential test \code{\link{stbp_composite}}.
#'
#'
#' @param data For count data, a numeric vector with for a single sampling bout
#' (NAs allowed). For binomial data, a matrix with observations in col 1 and
#' samples in col 2 (NAs \emph{not} allowed).
#' @param greater_than logical; if TRUE (default), the tested hypothesis is of
#' the form \eqn{H:\mu > \psi} otherwise, \eqn{H:\mu < \psi}.
#' @param hypothesis Single non-negative value with the hypothesized value
#' of \eqn{\mu}.
#' @param density_func Kernel probability density function for the data. See details.
#' @param overdispersion A character string (if a function) or a number
#' specifying the overdispersion parameter. Only required when using
#' \code{"negative binomial"} or \code{"beta-binomial"} as kernel densities.
#' See details.
#' @param prior Single number with initial prior. Must be on the interval \eqn{[0,1]}.
#' @param lower_bnd Single number indicating the lower bound of the parameter
#' space for \eqn{\mu}. Most cases is \eqn{0} (default).
#' @param upper_bnd Single number indicating the upper bound of the parameter
#' space for \eqn{\mu}. For count data, is often \code{Inf} (default), but it
#' must be \eqn{\leq 1} for binomial data.
#' @returns A single probability
#' @details
#' The \code{density_func} argument should be specified as character string.
#' Acceptable options are \code{"poisson"}, \code{"negative binomial"}, \code{"binomial"}
#' and \code{"beta-binomial"}. The overdispersion parameter for \code{"negative binomial"}
#' and \code{"beta-binomial"} can be either a constant or a function of the mean.
#' If a function, it should be specified as a character string with the name of
#' an existing function. For options of empirical functions to describe
#' overdispersion as a function of the mean see Binns et al. (2000).
#'
#' @references Binns, M.R., Nyrop, J.P. & Werf, W.v.d. (2000) \emph{Sampling and
#' monitoring in crop protection: the theoretical basis for developing practical
#' decision guides}. CABI Pub., Wallingford, Oxon, UK; New York, N.Y.
#'
#' Rincon, D.F., McCabe, I. & Crowder, D.W. (2025) Sequential testing of
#' complementary hypotheses about population density. \emph{Methods in Ecology
#' and Evolution}. <https://doi.org/10.1111/2041-210X.70053>
#'
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
#'                           density_func = "poisson",
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf) # returns 0.60630278
#'
#' # Same analysis but with a negative binomial kernel.
#' # Note that 'overdispersion' can either be a positive number or a function.
#'
#' stbp_posterior_composite(data = counts,
#'                           greater_than = TRUE,
#'                           hypothesis = 2,
#'                           density_func = "negative binomial",
#'                           overdispersion = 2,
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf) # returns 0.72558593
#' ## End (Not run)
#' @export

stbp_posterior_composite <- function(data,
                                     greater_than = TRUE,
                                     hypothesis,
                                     density_func,
                                     overdispersion = NA,
                                     prior,
                                     lower_bnd = 0,
                                     upper_bnd = Inf) {
  if(is.na(hypothesis) | hypothesis < 0)
    stop("ERROR: hypothesized values for mu should be non-negative and
         can't be NA.")

  if(is.na(match(density_func, c("poisson", "negative binomial", "binomial",
                                 "beta-binomial"))))
    stop("ERROR: only distributions poisson, negative binomial, binomial and
         beta-binomial are supported.")

  if(density_func == "negative binomial" | density_func == "beta-binomial") {
    if(is.na(overdispersion) | overdispersion < 0)
      stop("ERROR: a non-negative overdipersion parameter value must be
           provided.")
  }


  # data check

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }

  if(density_func == "binomial" | density_func == "beta-binomial") {

    if(!is.matrix(data) | ncol(data) != 2)
      stop("ERROR: binomial data should be presented as a two-column matrix,
           with observations in col 1 and samples in col 2.")

    if(any(is.na(data) | !is.wholenumber(data) | data < 0))
      stop("ERROR: binomial data must be non-negative integers and
           can't contain NAs.")

    if(hypothesis > 1)
      stop("ERROR: hypothesized values for binomial variables must be
           proportions.")

    if(upper_bnd > 1)
      stop("ERROR: possible values for mu should be within the interval [0,1].")
  }

  if(density_func == "poisson" | density_func == "negative binomial") {
    if(all(is.na(data)))
      stop("ERROR: sampling bouts should contain at least one datum.")

    if(any(!is.wholenumber(na.omit(data)) | na.omit(data) < 0))
      stop("ERROR: count data must be non-negative integers.")
  }

  # construction of likelihood function

  if(density_func == "poisson")
    likelihood_func <- function(data, x) dpois(data, lambda = x)


  if(density_func == "negative binomial") {
    if(is.character(overdispersion)) {
      likelihood_func <- function(data, x)
        dnbinom(data, mu = x,
                size = if(eval(parse(text =
                                     paste(overdispersion, "(x)"))) < 0 |
                          is.nan(eval(parse(text =
                                            paste(overdispersion, "(x)")))) |
                          is.na(eval(parse(text =
                                           paste(overdispersion, "(x)"))))) 0
                else eval(parse(text =
                                  paste(overdispersion, "(x)"))))
    }

    if(is.numeric(overdispersion)) {
      likelihood_func <- function(data, x)
        dnbinom(data, mu = x, size = overdispersion)
    }
  }

  if(density_func == "binomial")
    likelihood_func <- function(data, x)
      dbinom(data[,1], size = data[,2], prob = x)

  if(density_func == "beta-binomial") {
    if(is.character(overdispersion)) {
      likelihood_func <- function(data, x)
        emdbook::dbetabinom(data[,1], size = data[,2],
                            prob = x,
                            theta = if(eval(parse(text =
                                                  paste(overdispersion, "(x)"))) < 0 |
                                       is.nan(eval(parse(text =
                                                         paste(overdispersion, "(x)")))) |
                                       is.na(eval(parse(text =
                                                        paste(overdispersion, "(x)"))))) 0
                            else eval(parse(text =
                                              paste(overdispersion, "(x)"))))
    }

    if(is.numeric(overdispersion)) {
      likelihood_func <- function(data, x)
        emdbook::dbetabinom(data[,1],
                            size = data[,2],
                            prob = x,
                            theta = overdispersion)
    }
  }

  if(!is.matrix(data)) data <- as.numeric(na.omit(data))

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
#' @param data For count data, either a vector (for purely sequential designs)
#' or a matrix (group sequential designs) with sequential (non-negative) count
#' data, with sampling bouts collected over time in columns and samples within
#' bouts in rows. NAs are allowed in case sample size within bouts is unbalanced.
#' For binomial data, a list of matrices with integer non-negative values of
#' observations in col 1 and number of samples in col 2, so that each matrix
#' within the list corresponds to a sampling bout. NAs are \emph{not} allowed for
#' binomial data.
#' @param greater_than logical; if TRUE (default), the tested hypothesis is of
#' the form \eqn{H:\mu > \psi} otherwise, \eqn{H:\mu < \psi}.
#' @param hypothesis Either a single non-negative value or a vector of
#' non-negative values with the hypothesized population densities, \eqn{\mu}.
#' If a vector, it should contain at least as many values as \code{ncol(data)}
#' for count data or as \code{length(data)} for binomial data.
#' @param density_func Kernel probability density function for the data. See
#' details.
#' @param overdispersion A character string (if a function) or a number
#' specifying the overdispersion parameter. Only required when using
#' \code{"negative binomial"} or \code{"beta-binomial"} as kernel densities.
#' See details.
#' @param prior Single number with initial prior. Must be on the interval
#' \eqn{[0,1]}.
#' @param lower_bnd Single number indicating the lower bound of the parameter
#' space for \eqn{\mu}. Most cases is \eqn{0} (default).
#' @param upper_bnd Single number indicating the upper bound of the parameter
#' space for \eqn{\mu}. For count data, is often \code{Inf} (default), but it
#' must be \eqn{\leq 1} for binomial data.
#' @param lower_criterion Criterion to decide against the tested hypothesis.
#' This is the lowest credibility to the hypothesis to stop sampling and
#' decide against.
#' @param upper_criterion Criterion to decide in favor of the tested hypothesis.
#' This is the greatest credibility to the hypothesis to stop sampling and
#' decide in favor.
#'
#' @details
#' The \code{density_func} argument should be specified as character string.
#' Acceptable options are \code{"poisson"}, \code{"negative binomial"},
#' \code{"binomial"} and \code{"beta-binomial"}. The overdispersion
#' parameter for \code{"negative binomial"} and \code{"beta-binomial"} can be
#' either a constant or a function of the mean. If a function, it should be
#' specified as a character string with the name of an existing function.
#' For options of empirical functions to describe overdispersion as a function
#' of the mean see Binns et al. (2000).
#'
#' @references Binns, M.R., Nyrop, J.P. & Werf, W.v.d. (2000) \emph{Sampling and
#' monitoring in crop protection: the theoretical basis for developing practical
#' decision guides}. CABI Pub., Wallingford, Oxon, UK; New York, N.Y.
#'
#' Rincon, D.F., McCabe, I. & Crowder, D.W. (2025) Sequential testing of
#' complementary hypotheses about population density. \emph{Methods in Ecology
#' and Evolution}. <https://doi.org/10.1111/2041-210X.70053>
#'
#' @returns
#' An object of class \code{"STBP"}.
#'
#' @examples
#' # Testing the hypothesis of a population size being greater than 5 individuals
#' # per sampling unit (H: mu > 5). The sequential sampling is made of 5 sampling
#' # bouts made of one sample each.
#'
#' set.seed(101)
#' counts3 <- rpois(5, lambda = 3)
#'
#' test1F <- stbp_composite(data = counts3,
#'                           greater_than = TRUE,
#'                           hypothesis = 5,
#'                           density_func = "poisson",
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf,
#'                           lower_criterion = 0.001,
#'                           upper_criterion = 0.999)
#' test1F
#'
#' # returns "reject H".
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
#' test2F <- stbp_composite(data = countP,
#'                           greater_than = TRUE,
#'                           hypothesis = H,
#'                           density_func = "poisson",
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf,
#'                           lower_criterion = 0.001,
#'                           upper_criterion = 0.999)
#' test2F
#'
#' # returns "reject H".
#'
#' # Testing the hypothesis of a proportion of infested units being greater than
#' # 20% per sampling unit (H: mu > 0.2). The sequential sampling is made of 7
#' # sampling bouts each with 5 clusters of 10 samples from which binomial
#' # observations are recorded.
#'
#' set.seed(101)
#'
#' # binomial data generated with mu (prob) 0.05 over the hypothesized
#' # value (0.2)
#'
#' counts4 <- list()
#' for(i in 1: 7) {
#'   counts4[[i]] <- matrix(c(rbinom(5, size = 10, prob = 0.25), rep(10, 5)),
#'                         5, 2)
#' }
#'
#' # Run the test. Notice that upper_bnd = 1!
#'
#' test3F <- stbp_composite(data = counts4,
#'                           greater_than = TRUE,
#'                           hypothesis = 0.2,
#'                           density_func = "binomial",
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = 1,
#'                           lower_criterion = 0.001,
#'                           upper_criterion = 0.999)
#'
#' test3F # returns accept H after 3 sampling bouts
#'
#' ## End (Not run)
#' @export
#'
stbp_composite <- function(data,
                           greater_than = TRUE,
                           hypothesis,
                           density_func,
                           overdispersion = NA,
                           prior = 0.5,
                           lower_bnd = 0,
                           upper_bnd = Inf,
                           lower_criterion = 0.01,
                           upper_criterion = 0.99) {

  call <- match.call()

  if(density_func == "poisson" | density_func == "negative binomial") {

    # If hypothesis is just a single repeated value,
    # make a vector of that value repeated as many times as there are bouts.
    # This makes it so that the user can input either a hypothesis made of
    # a single value or a trajectory.

    if(is.vector(data))
      data <- matrix(data, 1, length(data))

    # Init vector with length equal to number of sampling bouts
    # and with initial prior as its first value

    if(length(hypothesis) == 1)
      hypothesis <- rep(hypothesis, ncol(data))

    if(length(hypothesis) < ncol(data))
      stop("ERROR: a hypothesized value for mu should be provided for each
           sampling bout.")

    posteriors <- c(prior, rep(NA, ncol(data) - 1))
    for(i in 1: ncol(data)) {
      if(posteriors[i] < 0.001) posteriors[i] <- 0.001
      if(posteriors[i] > 0.999) posteriors[i] <- 0.999
      bout = data[, i]
      posteriors[i + 1] = stbp_posterior_composite(bout,
                                                   greater_than,
                                                   hypothesis[i],
                                                   density_func,
                                                   overdispersion,
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
  }

  if(density_func == "binomial" | density_func == "beta-binomial") {

    if(!is.list(data))
      stop("ERROR: binomial sequential data should be provided as a list of
           matrices. See ?stbp_composite.")

    if(length(hypothesis) == 1)
      hypothesis <- rep(hypothesis, length(data))

    if(length(hypothesis) < length(data))
      stop("ERROR: a hypothesized value for mu should be provided for each
           sampling bout.")

    posteriors <- c(prior, rep(NA, length(data) - 1))
    for(i in 1: length(data)) {
      if(posteriors[i] < 0.001) posteriors[i] <- 0.001
      if(posteriors[i] > 0.999) posteriors[i] <- 0.999
      bout = data[[i]]
      posteriors[i + 1] = stbp_posterior_composite(bout,
                                                   greater_than,
                                                   hypothesis[i],
                                                   density_func,
                                                   overdispersion,
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
  }

  if(((posteriors[i + 1] < lower_criterion) ||
      (posteriors[i + 1] > upper_criterion))
  ) {
    if(posteriors[i + 1] < lower_criterion) response <- "reject H"
    else response <- "accept H"
  }
  else response <- "keep sampling"

  posteriors <- as.numeric(na.omit(posteriors))

  resp <- new("STBP",
              call = call,
              probabilities = posteriors,
              recommendation = response,
              iterations = i)

  resp
}


#' Posterior calculation for simple hypotheses about species absence
#'
#' @description
#' This function calculates a posterior probability for hypotheses about population
#' densities, of the form \eqn{H:\mu = 0}, given the data at a single
#' iteration. This function is to be used in a sequential framework, and called
#' on the sequential test \code{\link{stbp_simple}}.
#'
#' @param data For count data, a numeric vector with for a single sampling bout
#' (NAs allowed). For binomial data, a matrix with observations in col 1 and
#' samples in col 2 (NAs \emph{not} allowed).
#' @param density_func Kernel probability density function for the data. See details.
#' @param overdispersion A character string (if a function) or a number
#' specifying the overdispersion parameter. Only required when using
#' \code{"negative binomial"} or \code{"beta-binomial"} as kernel densities.
#' See details.
#' @param prior Single number with initial prior. Must be in the interval \eqn{[0,1]}.
#' @param upper_bnd Single number indicating the greatest possible value for \eqn{\mu}.
#' For count data, is often \code{Inf} (default), but it must be \eqn{\leq 1} for
#' binomial data.
#'
#' @return A single probability
#'
#' @details
#' The \code{density_func} argument should be specified as character string.
#' Acceptable options are \code{"poisson"}, \code{"negative binomial"},
#' \code{"binomial"} and \code{"beta-binomial"}. The overdispersion
#' parameter for \code{"negative binomial"} and \code{"beta-binomial"} can be
#' either a constant or a function of the mean. If a function, it should be
#' specified as a character string with the name of an existing function. For
#' options of empirical functions to describe overdispersion as a function of
#' the mean see Binns et al. (2000).
#'
#' @references Binns, M.R., Nyrop, J.P. & Werf, W.v.d. (2000) \emph{Sampling and
#' monitoring in crop protection: the theoretical basis for developing practical
#' decision guides}. CABI Pub., Wallingford, Oxon, UK; New York, N.Y.
#'
#' Rincon, D.F., McCabe, I. & Crowder, D.W. (2025) Sequential testing of
#' complementary hypotheses about population density. \emph{Methods in Ecology
#' and Evolution}. <https://doi.org/10.1111/2041-210X.70053>
#'
#' @examples
#' # Counts collected in a single sampling bout
#' counts <- c(0, 0, 0)
#'
#' # Calculate posterior probability from a naive 0.5 prior for H:mu=0
#' # (a species being absent in an area) with a poisson kernel.
#'
#' stbp_posterior_simple(data = counts,
#'                        density_func = "poisson",
#'                        prior = 0.5,
#'                        upper_bnd = Inf) # returns 0.75
#' ## End (Not run)
#' @export

stbp_posterior_simple <- function(data,
                                  density_func,
                                  overdispersion = NA,
                                  prior,
                                  upper_bnd = Inf) {

  if(is.na(match(density_func, c("poisson", "negative binomial", "binomial",
                                 "beta-binomial"))))
    stop("ERROR: only distributions poisson, negative binomial, binomial and
         beta-binomial are supported.")

  if(density_func == "negative binomial" | density_func == "beta-binomial") {
    if(is.na(overdispersion) | overdispersion < 0)
      stop("ERROR: a non-negative overdipersion parameter value must be
           provided.")
  }


  # data check

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }

  if(density_func == "binomial" | density_func == "beta-binomial") {

    if(!is.matrix(data) | ncol(data) != 2)
      stop("ERROR: binomial data should be presented as a two-column matrix,
           with observations in col 1 and samples in col 2.")

    if(any(is.na(data) | !is.wholenumber(data) | data < 0))
      stop("ERROR: binomial data must be non-negative integers and
           can't contain NAs.")

    if(upper_bnd > 1)
      stop("ERROR: possible values for mu should be within the interval [0,1].")
  }

  if(density_func == "poisson" | density_func == "negative binomial") {
    if(all(is.na(data)))
      stop("ERROR: sampling bouts should contain at least one datum.")

    if(any(!is.wholenumber(na.omit(data)) | na.omit(data) < 0))
      stop("ERROR: count data must be non-negative integers.")
  }

  # construction of likelihood functions

  if(density_func == "poisson")
    likelihood_func <- function(data, x) dpois(data, lambda = x)


  if(density_func == "negative binomial") {
    if(is.character(overdispersion)) {
      likelihood_func <- function(data, x)
        dnbinom(data, mu = x,
                size = if(eval(parse(text =
                                     paste(overdispersion, "(x)"))) < 0 |
                          is.nan(eval(parse(text =
                                            paste(overdispersion, "(x)")))) |
                          is.na(eval(parse(text =
                                           paste(overdispersion, "(x)"))))) 0
                else eval(parse(text =
                                  paste(overdispersion, "(x)"))))
    }

    if(is.numeric(overdispersion)) {
      likelihood_func <- function(data, x)
        dnbinom(data, mu = x, size = overdispersion)
    }
  }

  if(density_func == "binomial")
    likelihood_func <- function(data, x)
      dbinom(data[,1], size = data[,2], prob = x)

  if(density_func == "beta-binomial") {
    if(is.character(overdispersion)) {
      likelihood_func <- function(data, x)
        emdbook::dbetabinom(data[,1], size = data[,2],
                            prob = x,
                            theta = if(eval(parse(text =
                                                  paste(overdispersion, "(x)"))) < 0 |
                                       is.nan(eval(parse(text =
                                                         paste(overdispersion, "(x)")))) |
                                       is.na(eval(parse(text =
                                                        paste(overdispersion, "(x)"))))) 0
                            else eval(parse(text =
                                              paste(overdispersion, "(x)"))))
    }

    if(is.numeric(overdispersion)) {
      likelihood_func <- function(data, x)
        emdbook::dbetabinom(data[,1],
                            size = data[,2],
                            prob = x,
                            theta = overdispersion)
    }
  }

  if(!is.matrix(data)) data <- as.numeric(na.omit(data))

  likelihood <- function(x) {
    prod(likelihood_func(data, x))
  }

  H1 <- prior *
    likelihood(0)
  H0 <- (1 - prior) *
    integrate(
      Vectorize(likelihood),
      lower = 0,
      upper = Inf
    )$value
  posterior <- H1 / (H0 + H1)
  posterior
}

#' Sequential test of Bayesian posterior probabilities for simple hypotheses
#' about species absence
#'
#' @description
#' Runs a Sequential test of Bayesian Posterior Probabilities for hypotheses
#' about species absence of the form \eqn{H:\mu = 0}. Data is treated in a
#' sequential framework.
#'
#'
#' @param data For count data, either a vector (for purely sequential designs) o
#' a matrix (group sequential designs) with sequential count data, with sampling
#' bouts collected over time in columns and samples within bouts in rows. NAs
#' are allowed in case sample size within bouts is unbalanced. For binomial
#' data, a list of matrices with integer non-negative values of observations
#' in col 1 and number of samples in col 2, so that each matrix within the list
#' corresponds to a sampling bout. NAs are \emph{not} allowed for binomial data.
#' @param density_func Kernel probability density function for the data. See details.
#' @param overdispersion A character string (if a function) or a number
#' specifying the overdispersion parameter. Only required when using
#' \code{"negative binomial"} or \code{"beta-binomial"} as kernel densities.
#' See details.
#' @param prior Single number with initial prior. Must be in the interval \eqn{[0,1]}.
#' @param upper_bnd Single number indicating the greatest possible value for \eqn{\mu}.
#' For count data, is often \code{Inf} (default), but it must be \eqn{\leq 1}
#' for binomial data.
#' @param lower_criterion Criterion to decide against the tested hypothesis.
#' This is the lowest credibility to the hypothesis to stop sampling and decide
#' against.
#' @param upper_criterion Criterion to decide in favor of the tested hypothesis.
#' This is the greatest credibility to the hypothesis to stop sampling and decide
#' in favor.
#'
#' @details
#' The \code{density_func} argument should be specified as character string.
#' Acceptable options are \code{"poisson"}, \code{"negative binomial"},
#' \code{"binomial"} and \code{"beta-binomial"}. The overdispersion
#' parameter for \code{"negative binomial"} and \code{"beta-binomial"} can be
#' either a constant or a function of the mean. If a function, it should be
#' specified as a character string with the name of an existing function. For
#' options of empirical functions to describe overdispersion as a function of
#' the mean see Binns et al. (2000).
#'
#' @references Binns, M.R., Nyrop, J.P. & Werf, W.v.d. (2000) \emph{Sampling and
#' monitoring in crop protection: the theoretical basis for developing practical
#' decision guides}. CABI Pub., Wallingford, Oxon, UK; New York, N.Y.
#'
#' Rincon, D.F., McCabe, I. & Crowder, D.W. (2025) Sequential testing of
#' complementary hypotheses about population density. \emph{Methods in Ecology
#' and Evolution}. <https://doi.org/10.1111/2041-210X.70053>
#'
#' @returns
#' An object of class \code{"STBP"}.
#'
#' @examples
#' # Testing the absence of a species in a given area from a sequential random
#' # sampling of 3 bouts made of 10 samples (counts) each (all absences). Upper
#' # criterion set to 0.9999
#'
#' counts10 <- matrix(rep(0, 30), 10, 3)
#'
#' test1G <- stbp_simple(data = counts10,
#'                         density_func = "poisson",
#'                         prior = 0.5,
#'                         upper_bnd = Inf,
#'                         lower_criterion = 0,
#'                         upper_criterion = 0.9999)
#' test1G
#'
#' # returns a recommendation of "keep sampling" due to insufficient evidence.
#'
#' # Testing the same hypothesis with the same upper criterion but from a
#' # sequential random sampling of 3 bouts made of 30 samples (counts) each
#' # (all absences).
#'
#' counts30 <- matrix(rep(0, 90), 30, 3)
#'
#' test2G <- stbp_simple(data = counts30,
#'                         density_func= "poisson",
#'                         prior = 0.5,
#'                         upper_bnd = Inf,
#'                         lower_criterion = 0,
#'                         upper_criterion = 0.9999)
#' test2G
#'
#' # returns a recommendation of "accept H" of the species being absent from
#' # that area.
#'
#' ## End (Not run)
#' @export
stbp_simple <- function(data,
                        density_func,
                        overdispersion = NA,
                        prior = 0.5,
                        upper_bnd = Inf,
                        lower_criterion,
                        upper_criterion) {

  call <- match.call()

  if(density_func == "poisson" | density_func == "negative binomial") {

    if(is.vector(data))
      data <- matrix(data, 1, length(data))


    posteriors <- c(prior, rep(NA, ncol(data) - 1))
    for(i in 1: ncol(data)) {
      bout = data[, i]
      posteriors[i + 1] = stbp_posterior_simple(bout,
                                                density_func,
                                                overdispersion,
                                                prior = posteriors[i],
                                                upper_bnd = upper_bnd)
      # Break from iteration early if early_return is true,
      # the minimum iterations have been reached, and
      # if either of the decision criteria have been reached
      if(((posteriors[i + 1] < lower_criterion) ||
          (posteriors[i + 1] > upper_criterion))
      ) break
    }
  }

  if(density_func == "binomial" | density_func == "beta-binomial") {

    if(!is.list(data))
      stop("ERROR: binomial sequential data should be provided as a list of
           matrices. See ?stbp_simple.")

    posteriors <- c(prior, rep(NA, length(data) - 1))
    for(i in 1: length(data)) {
      bout = data[[i]]
      posteriors[i + 1] = stbp_posterior_simple(bout,
                                                density_func,
                                                overdispersion,
                                                prior = posteriors[i],
                                                upper_bnd = upper_bnd)
      # Break from iteration early if early_return is true,
      # the minimum iterations have been reached, and
      # if either of the decision criteria have been reached
      if(((posteriors[i + 1] < lower_criterion) ||
          (posteriors[i + 1] > upper_criterion))
      ) break
    }
  }



  if(((posteriors[i + 1] < lower_criterion) ||
      (posteriors[i + 1] > upper_criterion))
  ) {
    if(posteriors[i + 1] < lower_criterion) response <- "reject H"
    else response <- "accept H"
  }
  else response <- "keep sampling"

  posteriors <- as.numeric(na.omit(posteriors))

  resp <- new("STBP",
              call = call,
              probabilities = posteriors,
              recommendation = response,
              iterations = i)

  resp
}
