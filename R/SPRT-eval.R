#' Simulation of an SPRT
#'
#' @description
#' This function is called on Evaluation of Sequential probability ratio test,
#' \code{\link{SPRT.eval}}.
#'
#' @param pop_mean True population density.
#' @param obj A SPRT object
#' @param overdispersion.sim Overdispersion parameter or function for simulations.
#' @param seed Optional seed for random count generation.
#'
#' @returns
#' A list with average number of bouts required to reach a decision and
#' recommendation for H1
int_eval2 <- function(pop_mean, obj, overdispersion.sim = NA, seed = NULL){

  if (!is.null(seed)) set.seed(seed)
  samples <- rep(NA, 100)

  if(as.character(obj@call[5]) == "negative binomial") {

    if(is.na(overdispersion.sim)) {
      if(suppressWarnings(is.na(as.numeric(as.character(obj@call[6]))))) {
        overdispersion.sim <- eval(parse(text = as.character(obj@call[6])))
      } else
        overdispersion.sim <- as.numeric(as.character(obj@call[6]))
    }

    if(is.character(overdispersion.sim)) {
      pool <- rnbinom(n = 6000, mu = pop_mean,
                      size = if(eval(parse(text =
                                           paste0(overdispersion.sim, "(pop_mean)"))) < 0 |
                                is.nan(eval(parse(text = paste0(overdispersion.sim, "(pop_mean)")))) |
                                is.na(eval(parse(text = paste0(overdispersion.sim, "(pop_mean)"))))) 0
                      else eval(parse(text =
                                        paste0(overdispersion.sim, "(pop_mean)"))))
    }

    if(is.numeric(overdispersion.sim)) {
      pool <- rnbinom(mu = pop_mean, size = overdispersion.sim, n = 6000)
    }
  }

  if(as.character(obj@call[5]) == "poisson"){
    pool <- rpois(lambda = pop_mean, n = 6000)
  }

  if(as.character(obj@call[5]) == "binomial"){
    pool <- rbinom(prob = pop_mean, size = 1, n = 6000)
  }


  for(i in 1:100){
    samples[i] <- sample(pool, size = 1, replace = FALSE)
    pool <- pool[-match(samples[i], pool)]
  }


  a <- deparse(obj@call, width.cutoff = 20)
  a[1] <- gsub(as.character(obj@call[2]), "samples",
               deparse(obj@call, width.cutoff = 20)[1])


  ntest <- parse(text = a)
  test <- eval(ntest)

  return(list(
    samples = test@iterations,
    recommendation = test@recommendation
  ))
}


#' Evaluation of Sequential probability ratio tests
#'
#' @description
#' Obtains the average number of sampling bouts and the rate of acceptance for
#' \eqn{H_{1}}, from a sequential test of \eqn{H_{0}:\mu = \psi_{0}} vs.
#' \eqn{H_{1}:\mu < \psi_{1}} across a range of true population densities,
#' based on simulations. Sometimes called "operating characteristics".
#'
#' @param obj An object of class \code{"SPRT"}.
#' @param eval.range A vector with a sequence of true population densities to
#' evaluate.
#' @param overdispersion.sim A character string (if a function) or a non-negative number
#' specifying the overdispersion parameter used to generate simulated counts.
#' Only required when using \code{"negative binomial"} or \code{"beta-binomial"}
#' as kernel densities. See details.
#' @param N Number of simulations per true population density being evaluated.
#' @param seed Optional seed for random count generation.
#'
#' @details
#' The kernel probability density function to evaluate the test is that specified in the
#' argument \code{density_func} to create the \code{"SPRT"} object, but overdispersion
#' can be different to generate simulated counts. If \code{"negative binomial"} is used
#' as kernel density for the test and \code{overdispersion.sim} is not specified (NA),
#' then the same specification of the test is used to generate the counts. Ideally,
#' overdispersion for simulations should include uncertainty about the parameter to
#' produce more robust test evaluations. For example, if using a \code{negative binomial} kernel
#' and the Taylor's Power Law approach to obtain overdispersion, then overdispersion for
#' simulations should be specified as: \deqn{k = \frac{\mu^2}{a \mu^b e^z - \mu}} where \eqn{k}
#' is the overdispersion parameter of the negative binomial distribution, \eqn{a} and \eqn{b}
#' are parameters of the Taylor's Power Law and \eqn{z} is a normally distributed
#' variable with mean \eqn{0} and standard deviation \eqn{\sigma_{e}}, which is the root
#' of the mean square error for the regression used to estimate \eqn{a} and \eqn{b}. See
#' examples.
#'
#' @references Binns, M.R., Nyrop, J.P. & Werf, W.v.d. (2000) \emph{Sampling and
#' monitoring in crop protection: the theoretical basis for developing practical
#' decision guides}. CABI Pub., Wallingford, Oxon, UK; New York, N.Y.
#'
#' @returns
#' A list with the average number of sampling bouts required to reach a decision
#' (\code{$AvgSamples}), and the rate of acceptance for \eqn{H_{1}} across
#' the provided range of population densities (\code{$AcceptRate}).
#'
#' @export
#'
#' @examples
#'
#' # Function that describes negative binomial overdisperion from the mean
#' # and Taylor's Power Law parameters, a and b:
#'
#' estimate_k <- function(mean) {
#'                         a = 1.830012
#'                         b = 1.218041 # a and b are Taylor's Power Law parameters
#'                         (mean^2) / ((a * mean^(b)) - mean)
#'                        }
#'
#' # Generate a SPRT object with negative binomial and varying overdispersion
#'
#' counts <- c(2, 5, 6, 2, 7)
#'
#' test11 <- sprt(data = counts,
#'                 mu0 = 8,
#'                 mu1 = 10,
#'                 density_func = "negative binomial",
#'                 overdispersion = estimate_k(9),
#'                 alpha = 0.1,
#'                 beta = 0.1)
#'
#' # Stochastic version of 'estimate k', where sd here is the the root of the
#' # mean square error for the regression used to estimate a and b.
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
#' # Run model evaluation for test11 with varying overdispersion and
#' # added stochasticity.
#'
#' eval4 <- SPRT.eval(test11, eval.range = seq(4, 13),
#'                     overdispersion.sim = "estimate_k_stoch", N = 30)
#'
#' plot(seq(4, 13), eval4$AvgSamples, type = "o", xlab = "True population size",
#'       ylab = "Average number of bouts")
#' plot(seq(4, 13), eval4$AcceptRate, type = "o", xlab = "True population size",
#'       ylab = "Acceptance rate of H1")
#'
#' ## End (Not run)
SPRT.eval <- function(obj,
                      eval.range,
                      overdispersion.sim = NA,
                      N,
                      seed = NULL) {

  res <- list()
  for(i in 1:length(eval.range)) {
    res[[i]] <- replicate(N, int_eval2(pop_mean = eval.range[i], obj = obj,
                                       overdispersion.sim = overdispersion.sim,
                                       seed = seed))
  }

  av.samples <- rep(NA, length(eval.range))
  prop.acc <- rep(NA, length(eval.range))

  for(i in 1:length(eval.range)) {
    av.samples[i] <- mean(unlist(res[[i]][1, ]))
    prop.acc[i] <- length(which(unlist(res[[i]][2, ]) == "accept H1")) / N
  }

  return(list(AvgSamples = av.samples,
              AcceptRate = prop.acc))
}
