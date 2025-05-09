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
    if(suppressWarnings(is.na(as.numeric(as.character(object@call[4])))))
      psi <- paste("Trajectory = ", object@call[4]) else
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


#' Draws a plot with the sequence of posterior probabilities from a sequential
#' test of Bayesian posterior probabilities
#'
#' @description
#' Method for signature \code{"STBP"} to display resulting probabilities.
#'
#' @param x Created as a result of a call to \code{\link{stbp_simple}} or
#' \code{\link{stbp_composite}}.
#' @param y Unused entry
#'
#' @returns
#' A plot with the sequence of posterior probabilities
#'
#' @export
#'
#' @docType methods
#'
#' @examples
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
#' plot(test2F)
#'
#' ## End (Not run)
#' @aliases plot
setMethod("plot", signature = c(x = "STBP", y = "missing"), function(x, y = NULL) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar), add = TRUE)

  par(mar = c(5, 7, 2, 2))
  plot(x@probabilities, ylim = c(0, 1), type = "o", cex.lab = 2, yaxt = "n",
       xaxt = "n", xlab = "Sampling bout", ylab = "", cex.axis = 2, lwd = 2)
  axis(2, at = seq(0, 1, 0.2), cex.axis = 2, las = 2)
  axis(1, at = seq(1, (x@iterations + 1)), labels = seq(0, x@iterations),
       cex.axis = 2)
  abline(h = as.numeric(as.character(x@call[9])), lty = 2)
  abline(h = as.numeric(as.character(x@call[10])), lty = 2)
  title(ylab = "Posterior probability", cex.lab = 2, line = 4.5)
})
