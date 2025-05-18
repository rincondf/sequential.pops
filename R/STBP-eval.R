#' Simulation of an STBP
#'
#' @description
#' This function is called on Evaluation of Sequential tests of
#' Bayesian posterior probabilities, \code{\link{STBP.eval}}.
#'
#' @param pop_mean True population density.
#' @param prior Initial prior.
#' @param n Sample size within bouts.
#' @param obj A STBP object
#' @param overdispersion.sim Overdispersion parameter or function for simulations.
#' @param seed Optional seed for random count generation.
#'
#' @returns
#' A list with average number of bouts required to reach a decision and recommendation for H
int_eval <- function(pop_mean, prior, n, obj, overdispersion.sim = NA, seed = NULL){

  if(length(eval(parse(text = as.character(obj@call[4])))) > 1)
    pop_mean <- eval(parse(text = as.character(obj@call[4]))) * pop_mean


  if (!is.null(seed)) set.seed(seed)
  samples <- matrix(NA, n, 100)

  if(as.character(obj@call[5]) == "negative binomial" |
     as.character(obj@call[5]) == "beta-binomial") {

    if(is.na(overdispersion.sim)) {
      if(suppressWarnings(is.na(as.numeric(as.character(obj@call[6]))))) {
        overdispersion.sim <- as.character(obj@call[6])
      } else
        overdispersion.sim <- as.numeric(as.character(obj@call[6]))
    }
  }

  pool <- matrix(NA, 6000, length(pop_mean))

  if(as.character(obj@call[5]) == "negative binomial") {

    for(i in 1:length(pop_mean)) {
      if(is.character(overdispersion.sim)) {
        pool[, i] <- rnbinom(n = 6000, mu = pop_mean[i],
                             size = if(eval(parse(text =
                                                  paste0(overdispersion.sim, "(pop_mean[i])"))) < 0 |
                                       is.nan(eval(parse(text = paste0(overdispersion.sim, "(pop_mean[i])")))) |
                                       is.na(eval(parse(text = paste0(overdispersion.sim, "(pop_mean[i])"))))) 0
                             else eval(parse(text =
                                               paste0(overdispersion.sim, "(pop_mean[i])"))))
      }

      if(is.numeric(overdispersion.sim)) {
        pool[, i] <- rnbinom(mu = pop_mean[i], size = overdispersion.sim, n = 6000)
      }
    }
  }

  if(as.character(obj@call[5]) == "beta-binomial"){

    for(i in 1:length(pop_mean)) {
      if(is.character(overdispersion.sim)) {
        pool[, i] <- rbetabinom(n = 6000, prob = pop_mean[i], size = 1,
                                theta = if(eval(parse(text =
                                                      paste0(overdispersion.sim, "(pop_mean[i])"))) < 0 |
                                           is.nan(eval(parse(text = paste0(overdispersion.sim, "(pop_mean[i])")))) |
                                           is.na(eval(parse(text = paste0(overdispersion.sim, "(pop_mean[i])"))))) 0
                                else eval(parse(text =
                                                  paste0(overdispersion.sim, "(pop_mean[i])"))))
      }

      if(is.numeric(overdispersion.sim)) {
        pool[, i] <- rbetabinom(prob = pop_mean[i], size = 1, theta = overdispersion.sim, n = 6000)
      }
    }
  }

  if(as.character(obj@call[5]) == "poisson"){

    for(i in 1:length(pop_mean)){
      pool[, i] <- rpois(lambda = pop_mean[i], n = 6000)
    }
  }

  if(as.character(obj@call[5]) == "binomial"){

    for(i in 1:length(pop_mean)){
      pool[, i] <- rbinom(prob = pop_mean[i], size = 1, n = 6000)
    }
  }


  if(length(pop_mean) == 1) {
    if(as.character(obj@call[5]) == "poisson" | as.character(obj@call[5]) == "negative binomial") {
      samples <- matrix(NA, n, 100)
      for(i in 1:100){
        samples[, i] <- sample(pool[!is.na(pool[, 1]), 1], size = n, replace = FALSE)
        pool[match(samples[, i], pool[, 1]), 1] <- NA
      }
    }


    if(as.character(obj@call[5]) == "binomial" | as.character(obj@call[5]) == "beta-binomial") {
      samples <- list()
      for(i in 1:100){
        samples[[i]] <- matrix(c(sample(pool[!is.na(pool[, 1]), 1], size = n, replace = FALSE), rep(1, n)), n, 2)
        pool[match(samples[[i]][, 1], pool[, 1]), 1] <- NA
      }
    }
  }

  if(length(pop_mean) > 1) {

    if(as.character(obj@call[5]) == "poisson" | as.character(obj@call[5]) == "negative binomial") {

      samples <- matrix(NA, n, length(pop_mean))
      for(i in 1:length(pop_mean)){
        samples[, i] <- sample(pool[, i], size = n, replace = FALSE)
      }
    }


    if(as.character(obj@call[5]) == "binomial" | as.character(obj@call[5]) == "beta-binomial") {
      samples <- list()
      for(i in 1:length(pop_mean)){
        samples[[i]] <- matrix(c(sample(pool[, i], size = n, replace = FALSE), rep(1, n)), n, 2)
      }
    }
  }


  a <- deparse(obj@call, width.cutoff = 20)
  a[1] <- "stbp_composite(data = samples, "
  ntest <- parse(text = a)

  test <- eval(ntest)

  return(list(
    samples = test@iterations,
    recommendation = test@recommendation
  ))
}

#' Evaluation of Sequential tests of Bayesian posterior probabilities
#'
#' @description
#' Obtains the average number of sampling bouts and the rate of acceptance for
#' \eqn{H:\mu > \psi} or \eqn{H:\mu < \psi} across a range of true population
#' densities, based on simulations. Sometimes called "operating characteristics".
#'
#' @param obj An object of class \code{"STBP"}.
#' @param eval.range The evaluation range. A vector with a sequence of true population
#' densities to evaluate. If the \code{hypothesis} argument in the test is a trajectory
#' (i.e., a vector), the evaluation range is a factor of the hypothesized trajectory.
#' See details.
#' @param n Sample size within bouts.
#' @param prior Single number with initial prior. Must be on the interval
#' \eqn{[0,1]}.
#' @param overdispersion.sim A character string (if a function) or a a non-negative number
#' specifying the overdispersion parameter used to generate simulated counts.
#' Only required when using \code{"negative binomial"} or \code{"beta-binomial"}
#' as kernel densities. See details.
#' @param N Number of simulations per true population density or trajectory being evaluated.
#' @param seed Optional seed for random count generation.
#'
#' @details
#' The kernel probability density function to evaluate the test is that specified in the
#' argument \code{density_func} to create the \code{"STBP"} object, but overdispersion
#' can be different to generate simulated counts. If \code{"negative binomial"} or
#' \code{"beta-binomial"} are used as kernel densities for the test and
#' \code{overdispersion.sim} is not specified (NA), then the same specification of the test
#' is used to generate the counts. Ideally, overdispersion for simulations should include
#' uncertainty about the parameter to produce more robust test evaluations. For example,
#' if using a \code{negative binomial} kernel and the Taylor's Power Law approach
#' to obtain overdispersion, then overdispersion for simulations should be specified
#' as: \deqn{k = \frac{\mu^2}{a \mu^b e^z - \mu}} where \eqn{k} is the overdispersion
#' parameter of the negative binomial distribution, \eqn{a} and \eqn{b} are parameters
#' of the Taylor's Power Law and \eqn{z} is a normally distributed variable with mean \eqn{0}
#' and standard deviation \eqn{\sigma_{e}}, which is the root of the mean square error
#' for the regression used to estimate \eqn{a} and \eqn{b}. See examples.
#'
#' The evaluation range in the \code{eval.range} argument should cover relevant population
#' densities for which the test will be applied. These densities are often around the
#' threshold to check sampling sizes and error rates at those critical levels. In the case
#' of dynamic hypotheses (vectors), the range should be given as factors of the hypothesized
#' trajectory. For example, to evaluate this test for a hypothesis trajectory the
#' \code{eval.range} argument could be \code{seq(0.1, 2, 0.1)}, so the evaluation runs from
#' \eqn{10} percent of the trajectory to twice the trajectory in \eqn{10} percent increments.
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
#' A list with the average number of sampling bouts required to reach a
#' decision (\code{$AvgSamples}), and the rate of acceptance for \eqn{H} across
#' the provided range of population densities (\code{$AcceptRate}).
#'
#' @export
#'
#' @examples
#'
#' # Assuming a negative binomial count variable whose overdispersion parameter,
#' # k, varies as a function of the mean, and that the variance-mean relationship
#' # is well described with Taylor's Power Law, a function to obtain k can be:
#'
#' estimate_k <- function(mean) {
#'                         a = 1.830012
#'                         b = 1.218041 # a and b are Taylor's Power Law parameters
#'                         (mean^2) / ((a * mean^(b)) - mean)
#'                         }
#'
#' # Generate some counts to create an STBP object with the model specifications
#'
#' counts3 <- rnbinom(20, mu = 5, size = estimate_k(5))
#'
#' # Run the test to create the STBP object
#'
#' test1F <- stbp_composite(data = counts3,
#'                           greater_than = TRUE,
#'                           hypothesis = 9,
#'                           density_func = "negative binomial",
#'                           overdispersion = "estimate_k",
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf,
#'                           lower_criterion = 0.01,
#'                           upper_criterion = 0.99)
#'
#' test1F
#'
#' # Model evaluation is carried out based on simulated counts, and more realistic
#' # counts could be generated if uncertainty about the overdispersion parameter is
#' # considered. A function to obtain values for the overdispersion parameter, k,
#' # with added stochasticity could be (following Binns et al. 2000):
#'
#' estimate_k_stoch <- function(mean) {
#'                               a <- 1.830012
#'                               b <- 1.218041
#'                               (mean^2) /
#'                                 ((a * mean^(b) *
#'                                   exp(truncdist::rtrunc(1,
#'                                                         "norm",
#'                                                         a = log(1 / (a * mean^(b - 1))),
#'                                                         b = Inf,
#'                                                         mean = 0,
#'                                                         sd = 0.3222354)))
#'                                          - mean)
#'                                    }
#'
#' # where sd here is the the root of the mean square error for the regression
#' # used to estimate a and b. Note that this is a stochastic version
#' # of 'estimate_k'.
#'
#' # Run model evaluation for testF1 with varying overdispersion and
#' # added stochasticity.
#'
#' eval1 <- STBP.eval(test1F,
#'                     eval.range = seq(2, 11),
#'                     n = 1, prior = 0.5,
#'                     overdispersion.sim = "estimate_k_stoch",
#'                     N = 20)
#'
#' plot(seq(2, 11), eval1$AvgSamples, type = "o", xlab = "True population size",
#'       ylab = "Average number of bouts")
#'
#' plot(seq(2, 11), eval1$AcceptRate, type = "o", xlab = "True population size",
#'       ylab = "Acceptance rate of H")
#'
#' # Alternatively, the evaluation could be carried out omitting variation about
#' # overdispersion. For that the overdispersion argument is omitted and the same
#' # specification of the model is used
#'
#' eval2 <- STBP.eval(test1F,
#'                     eval.range = seq(2, 11),
#'                     n = 1, prior = 0.5,
#'                     N = 20)
#'
#' plot(seq(2, 11), eval2$AvgSamples, type = "o", xlab = "True population size",
#'       ylab = "Average number of bouts")
#'
#' plot(seq(2, 11), eval2$AcceptRate, type = "o", xlab = "True population size",
#'       ylab = "Acceptance rate of H")
#'
#' # When there is no overdispersion (poisson or binomial distributions) the
#' # procedure is much simpler
#'
#' test2F <- stbp_composite(data = counts3,
#'                           greater_than = TRUE,
#'                           hypothesis = 5,
#'                           density_func = "poisson",
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf,
#'                           lower_criterion = 0.01,
#'                           upper_criterion = 0.99)
#'
#' test2F
#'
#' # Overdispersion is omitted here
#'
#' eval3 <- STBP.eval(test2F,
#'                   eval.range = seq(1, 8),
#'                   n = 1,
#'                   prior = 0.5, N = 20)
#'
#' plot(seq(1, 8), eval3$AvgSamples, type = "o", xlab = "True population size",
#'       ylab = "Average number of bouts")
#'
#' plot(seq(1, 8), eval3$AcceptRate, type = "o", xlab = "True population size",
#'       ylab = "Acceptance rate of H")
#'
#' # Variations if n, the sample size within each bout, can also be changed
#' # (not possible in SPRT)!
#'
#' ## End (Not run)
STBP.eval <- function(obj,
                      eval.range,
                      n,
                      prior,
                      overdispersion.sim = NA,
                      N,
                      seed = NULL) {

  res <- list()
  for(i in 1:length(eval.range)) {
    res[[i]] <- replicate(N, int_eval(pop_mean = eval.range[i], prior = prior,
                                      n = n, obj = obj, overdispersion.sim = overdispersion.sim,
                                      seed = seed))
  }

  av.samples <- rep(NA, length(eval.range))
  prop.acc <- rep(NA, length(eval.range))

  for(i in 1:length(eval.range)) {
    av.samples[i] <- mean(unlist(res[[i]][1, ]))
    prop.acc[i] <- length(which(unlist(res[[i]][2, ]) == "accept H")) / N
  }

  return(list(AvgSamples = av.samples,
              AcceptRate = prop.acc))
}
