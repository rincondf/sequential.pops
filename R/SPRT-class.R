#' Class "SPRT". Result of a Sequential Probability Ratio Test
#'
#' @description
#' This class encapsulates results of a Sequential Probability Ratio Test.
#'
#'
#' @slot call (language) The call to \code{\link{sprt}}.
#' @slot data (numeric) Vector with cumulative counts.
#' @slot hi_int (numeric) Intercept of the lower stop line.
#' @slot low_int (numeric) Intercept of the higher stop line.
#' @slot slope (numeric) Slope of stop lines.
#' @slot recommendation (character) Recommendation on whether to accept H0 or H1,
#' or keep sampling.
#' @slot iterations (numeric) Number of sequential sampling bouts required or processed.
#' @export
#'
#' @examples
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
#' show(test11) # returns "accept H1" after 5 sampling bouts processed.
#'
#' ## End (Not run)
setClass("SPRT",
         slots = c(
           call = "language",
           data = "numeric",
           hi_int = "numeric",
           low_int = "numeric",
           slope = "numeric",
           recommendation = "character",
           iterations = "numeric"
         ))
