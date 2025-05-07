#' Displays results from a "STBP" object nicely
#'
#' @description
#' Method for signature \code{"STBP"} to show results.
#'
#'
#' @param object Created as a result of a call to \code{\link{stbp_simple}} or \code{\link{stbp_composite}}.
#'
#' @returns
#' A summary of the test.
#'
#' @export
#'
#' @docType methods
#'
#' @examples
#' set.seed(101)
#' counts3 <- rpois(5, lambda = 3)
#'
#' test1F <- stbp_composite(data = counts3,
#'                           greater_than = TRUE,
#'                           hypothesis = 5,
#'                           likelihood_func = function(data, x)
#'                               dpois(data, lambda = x),
#'                           prior = 0.5,
#'                           lower_bnd = 0,
#'                           upper_bnd = Inf,
#'                           lower_criterion = 0.001,
#'                           upper_criterion = 0.999)
#' show(test1F)
#' # returns "reject H".
#'
#' counts10 <- matrix(rep(0, 30), 10, 3)
#'
#' test1G <- stbp_simple(data = counts10,
#'                         hypothesis = 0,
#'                         likelihood_func= function(data, x)
#'                             dpois(data, lambda = x),
#'                         prior = 0.5,
#'                         upper_bnd = Inf,
#'                         lower_criterion = 0,
#'                         upper_criterion = 0.9999)
#' show(test1G)
#'
#' # returns a recommendation of "keep sampling" due to insufficient evidence.
#'
#' @aliases show
setMethod("show", "STBP", function(object){

  dists <- c("dpois", "dnbinom", "dbinom", "dbetabinom")

  if(as.character(object@call[1]) == "stbp_simple") {
    sign <- "="
    psi <- as.numeric(as.character(object@call[3]))
    dist1 <- names(which(sapply(dists, grepl, as.character(object@call[4])) == TRUE))
    if(dist1 == "dpois") dist1 <- "Poisson"
    if(dist1 == "dnbinom") dist1 <- "Negative binomial"
    if(dist1 == "dbinom") dist1 <- "Binomial"
    if(dist1 == "dbetabinom") dist1 <- "Beta-binomial"
  }
  else {
    psi <- as.numeric(as.character(object@call[4]))
    if(as.character(object@call[3]) == TRUE) sign <- ">"
    if(as.character(object@call[3]) == FALSE) sign <- "<"
    dist1 <- names(which(sapply(dists, grepl, as.character(object@call[5])) == TRUE))
    if(dist1 == "dpois") dist1 <- "Poisson"
    if(dist1 == "dnbinom") dist1 <- "Negative binomial"
    if(dist1 == "dbinom") dist1 <- "Binomial"
    if(dist1 == "dbetabinom") dist1 <- "Beta-binomial"
  }

  cat("\nSequential test of Bayesian posterior probabilities\n")
  cat("Family:", dist1)
  cat("\nH: mu", sign, psi)
  cat("\nProbability: ")
  cat(round(as.numeric(object@probabilities[3]), 5),"from", object@iterations, "sampling bouts")
  cat("\nRecommendation based on provided criteria: ")
  cat(as.character(object@recommendation))
})
