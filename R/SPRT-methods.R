#' Displays results from a "SPRT" object nicely
#'
#' @description
#' Method for signature \code{"SPRT"} to show results.
#'
#'
#' @param object Created as a result of a call to \code{\link{sprt}}.
#'
#' @returns
#' A summary of the test.
#'
#' @export
#'
#' @docType methods
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
setMethod("show", "SPRT", function(object){
  cat("\nSequential Probability Ratio Test\n")
  cat("Family:", as.character(object@call[5]))
  cat("\nH0: mu =", as.character(object@call[3]))
  cat("\nH1: mu =", as.character(object@call[4]), "\n")
  cat("\nRecommendation based on provided criteria: ")
  cat(as.character(object@recommendation))
  cat("\nNumber of sampling bouts processed: ", object@iterations)
})



#' Draws a plot with the cumulative introduced counts contrasted with stop lines
#'
#' @description
#' Method for signature \code{"SPRT"} to display results.
#'
#'
#' @param x Created as a result of a call to \code{\link{sprt}}.
#' @param y Unused entry.
#'
#' @returns
#' A plot with cumulative counts contrasted with stop lines from a \code{"SPRT"}.
#'
#' @export
#'
#' @docType methods
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
#' plot(test11)
#'
#' ## End (Not run)
setMethod("plot", signature = c(x = "SPRT", y = "missing"), function(x, y) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar), add = TRUE)

  if(x@recommendation == "keep sampling") max.n <- length(x@data) + 5
  if(x@recommendation == "accept H0" | x@recommendation == "accept H1")
    max.n <- length(x@data) + 1

  par(mar = c(5, 7, 2, 2))
  plot(seq(1, length(x@data)), x@data, type = "o",
       ylim = c(0, max(x@data) + 1), xlim = c(0, max.n), cex.lab = 2, yaxt = "n",
       xaxt = "n", xlab = "Sampling bout",
       ylab = "", cex.axis = 2, lwd = 2, col = "blue")

  axis(2, at = seq(0, max(x@data) + 1, round((max(x@data) + 1)/5)), cex.axis = 2, las = 2)
  axis(1, at = seq(0, max.n, round(max.n/7)), cex.axis = 2)

  abline(a = x@hi_int, b = x@slope, lwd = 2)
  abline(a = x@low_int, b = x@slope, lwd = 2)

  title(ylab = "Cumulative counts", cex.lab = 2, line = 4.5)

})
