#' Sequential Probability Ratio Test
#'
#' @description
#' Runs a Sequential Probability Ratio Test for hypotheses
#' about population densities of the form \eqn{H_{0}:\mu = \psi_{0}} vs.
#' \eqn{H_{1}:\mu < \psi_{1}}, where \eqn{\psi_{0} < \psi_{1}}. Data is treated
#' in a sequential framework.
#'
#' @param data Optional vector of count data (NAs not allowed). Each value is
#' considered a sampling bout over time. Can't process group sequential data.
#' If not provided (NA), returns a \code{"SPRT"} object that can be used to plot
#' a chart with stop lines or to return stop lines coefficients. If provided, returns
#' a \code{"SPRT"} object that can be used to plot a chart with stop lines and data
#' or to return test summary.
#' @param mu0 Single non-negative number with the value for the low
#' hypothesized population density, \eqn{\psi_{0}}.
#' @param mu1 Single non-negative number with the value for the high
#' hypothesized population density, \eqn{\psi_{1}}.
#' @param density_func Kernel probability density function for the data. See
#' details.
#' @param overdispersion A number specifying the overdispersion parameter.
#' Only required when using \code{"negative binomial"} as kernel density.
#' See details.
#' @param alpha Single number indicating tolerable type I error rate.
#' @param beta Single number indicating tolerable type II error rate.
#'
#' @details
#' The \code{density_func} argument should be specified as character string.
#' Acceptable options are \code{"poisson"}, \code{"negative binomial"}, and
#' \code{"binomial"}. As far as we know, no one has ever calculated and published stop
#' lines for the beta-binomial family. The overdispersion parameter should only be
#' specified as a constant. In contrast to the STBP, SPRT is only use overdispersion
#' to calculate stop lines, so the estimate for the threshold population density
#' should be used (e.g., at \eqn{\psi_{0} < \psi < \psi_{1}}).
#'
#' @references Binns, M.R., Nyrop, J.P. & Werf, W.v.d. (2000) \emph{Sampling and
#' monitoring in crop protection: the theoretical basis for developing practical
#' decision guides}. CABI Pub., Wallingford, Oxon, UK; New York, N.Y.
#'
#' Wald, A. (1945) Sequential Tests of Statistical Hypotheses. \emph{The Annals
#' of Mathematical Statistics} 16(2): 117-186.
#'
#' @returns
#' An object of class \code{"SPRT"}.
#'
#' @export
#'
#' @examples
#'
#' # If no data is provided, an object of class "SPRT" is returned from which a
#' # chart with stop lines or a summary of the test with coefficients for stop lines
#' # can be extracted.
#'
#' test00 <- sprt(mu0 = 2,
#'                 mu1 = 4,
#'                 density_func = "negative binomial",
#'                 overdispersion = 4.6,
#'                 alpha = 0.1,
#'                 beta = 0.1)
#'
#' test00 # returns test specification and stop lines coefficients
#' plot(test00) # returns a chart with stop lines
#'
#'
#' # If data is provided, an object of class "SPRT" is returned with test results.
#'
#' counts <- c(2, 5, 6, 2, 7)
#'
#' test11 <- sprt(data = counts,
#'                 mu0 = 2,
#'                 mu1 = 4,
#'                 density_func = "negative binomial",
#'                 overdispersion = 4.6,
#'                 alpha = 0.1,
#'                 beta = 0.1)
#' test11 # returns "accept H1" after 5 sampling bouts processed.
#'
#' ## End (Not run)
#'
sprt <- function(data = NA,
                 mu0,
                 mu1,
                 density_func,
                 overdispersion,
                 alpha,
                 beta) {
  if(length(mu0) > 1 | length(mu1) > 1 | mu0 < 0 | mu1 < 0)
    stop("ERROR: hypothesized values, mu1 and mu0, should be single
         non-negative values.")

  if(is.na(match(density_func, c("poisson", "negative binomial", "binomial"))))
    stop("ERROR: only distributions poisson, negative binomial and binomial are
         supported.")

  if(density_func == "negative binomial") {
    if(is.na(overdispersion) | overdispersion < 0)
      stop("ERROR: a non-negative overdipersion parameter value must be
           provided.")
  }

  if(length(alpha) > 1 | length(beta) > 1 | alpha < 0 |
     beta < 0 | alpha > 1 | beta > 1)
    stop("ERROR: error rates should single values be within the interval [0,1].")


  call <- rlang::call_match(defaults = TRUE)

  if(density_func == "poisson") {

    low_int_pois <- function(al, be, me0, me1) {
      (log(be / (1 - al))) / (log((me1) / (me0)))
    }

    low_um <- low_int_pois(al = alpha, be = beta, me0 = mu0, me1 = mu1)

    hi_int_pois <- function(al, be, me0, me1) {
      (log((1 - al) / be)) / (log((me1) / (me0)))
    }

    hi_um <- hi_int_pois(al = alpha, be = beta, me0 = mu0, me1 = mu1)

    pen_pois <- function(me0, me1){
      (me1 - me0) / (log(me1 / me0))
    }

    pen <- pen_pois(me0 = mu0, me1 = mu1)
  }


  if(density_func == "binomial") {

    low_int_bin <- function(al, be, p1, p2) {
      q1 <- 1 - p1
      q2 <- 1 - p2
      (log(be / (1 - al))) / (log((p2*q1) / (p1*q2)))
    }

    low_um <- low_int_bin(al = alpha, be = beta, p1 = mu0, p2 = mu1)

    hi_int_bin <- function(al, be, p1, p2) {
      q1 <- 1 - p1
      q2 <- 1 - p2
      (log((1 - be) / (al))) / (log((p2*q1) / (p1*q2)))
    }

    hi_um <- hi_int_bin(al = alpha, be = beta, p1 = mu0, p2 = mu1)

    pen_bin <- function(p1, p2){
      q1 <- 1 - p1
      q2 <- 1 - p2
      (log(q1 / q2)) / (log((p2*q1) / (p1*q2)))
    }

    pen <- pen_bin(p1 = mu0, p2 = mu1)
  }


  if(density_func == "negative binomial") {

    low_int_nb <- function(al, be, me0, me1, k_est){
      (log(be / (1 - al))) / (log((me1 * (me0 + k_est)) / (me0 * (me1 + k_est))))
    }

    low_um <- low_int_nb(al = alpha, be = beta, me0 = mu0, me1 = mu1,
                         k_est = overdispersion)

    hi_int_nb <- function(al, be, me0, me1, k_est){
      (log((1 - be) / (al))) / (log((me1 * (me0 + k_est)) / (me0 * (me1 + k_est))))
    }

    hi_um <- hi_int_nb(al = alpha, be = beta, me0 = mu0, me1 = mu1,
                       k_est = overdispersion)

    ll_sl_nb <- function(al, be, me0, me1, k_est){
      (k_est * log((me1 + k_est) / (me0 + k_est))) /
        (log((me1 * (me0 + k_est)) / (me0 * (me1 + k_est))))
    }

    pen <- ll_sl_nb(al = alpha, be = beta, me0 = mu0, me1 = mu1,
                    k_est = overdispersion)
  }


  # Lower stop line:

  low_seq_c <- function(x){
    pen * x + low_um
  }

  # Upper stop line:

  hi_seq_c <- function(x){
    pen * x + hi_um
  }

  # if there no data, returns an object with test specification
  if(all(is.na(data))) {
    resp <- new("SPRT",
                call = call,
                data = NA,
                hi_int = hi_um,
                low_int = low_um,
                slope = pen,
                recommendation = NA,
                iterations = 0)
  }

  if(all(!is.na(data))) {

    if(!is.numeric(data))
      stop("ERROR: data should be numeric.")

    if(any(is.na(data)))
      stop("ERROR: no NAs allowed.")

    for(i in 1:length(data)){
      dat <- cumsum(data)[i]
      if((dat < low_seq_c(i)) | (dat > hi_seq_c(i))) break
    }

    if((dat < low_seq_c(i)) | (dat > hi_seq_c(i))) {
      if(dat < low_seq_c(i)) rec <- "accept H0"
      if(dat > hi_seq_c(i)) rec <- "accept H1"
    } else rec <- "keep sampling"

    resp <- new("SPRT",
                call = call,
                data = cumsum(data),
                hi_int = low_um,
                low_int = hi_um,
                slope = pen,
                recommendation = rec,
                iterations = i)
  }
  resp
}
