#' Class "STBP". Result of a Sequential tests of Bayesian posterior probabilities
#'
#' @description
#' This class encapsulates results of a Sequential tests of Bayesian posterior probabilities
#'
#' @slot call (language) The call to \code{\link{stbp_simple}} or \code{\link{stbp_composite}}.
#' @slot probabilities (numeric) Vector of sequential posterior probabilities.
#' @slot recommendation (character) Recommendation on H, whether to accept, reject or keep sampling.
#' @slot iterations (numeric) Number of sequential sampling bouts required or processed.
#'
#' @export
#'
#' @examples
#'
#' set.seed(101)
#' counts3 <- rpois(5, lambda = 3)
#'
#' test1F <- stbp_composite(data = counts3,
#'                           greater_than = TRUE,
#'                           hypothesis = 5,
#'                           likelihood_func = function(data, x)
#'                               {dpois(data, lambda = x)},
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf,
#'                           lower_criterion = 0.001,
#'                           upper_criterion = 0.999)
#' test1F
#' # returns "reject H".
#'
#' counts10 <- matrix(rep(0, 30), 10, 3)
#'
#' test1G <- stbp_simple(data = counts10,
#'                         hypothesis = 0,
#'                         likelihood_func= function(data, x)
#'                             {dpois(data, lambda = x)},
#'                         prior = 0.5,
#'                         upper_bnd = Inf,
#'                         lower_criterion = 0,
#'                         upper_criterion = 0.9999)
#' test1G
#'
#' # returns a recommendation of "keep sampling" due to insufficient evidence.
#'
setClass("STBP",
         slots = c(
           call = "language",
           probabilities = "numeric",
           recommendation = "character",
           iterations = "numeric"
         ))
