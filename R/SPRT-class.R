#' Class "SPRT". Result of a Sequential Probability Ratio Test
#'
#' @description
#' This class encapsulates results or specification of a Sequential Probability
#' Ratio Test.
#'
#'
#' @slot call (language) The call to \code{\link{sprt}}.
#' @slot data (numeric, logical) Either a vector with cumulative counts or NA.
#' @slot hi_int (numeric) Intercept of the lower stop line.
#' @slot low_int (numeric) Intercept of the higher stop line.
#' @slot slope (numeric) Slope of stop lines.
#' @slot recommendation (character, logical) Either a recommendation on whether
#' to accept H0 or H1, or keep sampling, or NA.
#' @slot iterations (numeric) Number of sequential sampling bouts required or
#' processed.
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
           data = "ANY",
           hi_int = "numeric",
           low_int = "numeric",
           slope = "numeric",
           recommendation = "ANY",
           iterations = "numeric"
         ))
