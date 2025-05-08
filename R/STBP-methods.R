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
#'                           density_func = "poisson",
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
#'                         density_func= "poisson",
#'                         prior = 0.5,
#'                         upper_bnd = Inf,
#'                         lower_criterion = 0,
#'                         upper_criterion = 0.9999)
#' show(test1G)
#'
#' # returns "keep sampling" due to insufficient evidence.
#'
#' ## End (Not run)
#' @aliases show
setMethod("show", "STBP", function(object){

  dists <- c("poisson", "negative binominal", "binomial", "beta-binomial")

  if(as.character(object@call[1]) == "stbp_simple") {
    sign <- "="
    psi <- 0
    dist1 <- as.character(object@call[4])
  }
  else {
    psi <- as.numeric(as.character(object@call[4]))
    if(as.character(object@call[3]) == TRUE) sign <- ">"
    if(as.character(object@call[3]) == FALSE) sign <- "<"
    dist1 <- as.character(object@call[5])
  }

  cat("\nSequential test of Bayesian posterior probabilities\n")
  cat("Family:", dist1)
  cat("\nH: mu", sign, psi)
  cat("\nProbability: ")
  cat(round(as.numeric(object@probabilities[length(object@probabilities)]), 5),
      "from", object@iterations, "sampling bouts")
  cat("\nRecommendation based on provided criteria: ")
  cat(as.character(object@recommendation))
})
