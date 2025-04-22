#' Set Optimized Plotting Parameters
#'
#' @description
#' Sets up optimized graphical parameters for R plots, particularly suited for
#' 3x3 inch square plots. This function provides a wrapper around base R's
#' par() with carefully chosen defaults for scientific visualization.
#'
#' @param mai Numeric vector of length 4 specifying margin size in inches
#'   (bottom, left, top, right). Default is c(0.9, 0.9, 0.9, 0.9).
#' @param las Numeric in {0,1,2,3}; style of axis labels. Default is 1
#'   (always horizontal).
#' @param tcl Numeric specifying length of tick marks as fraction of height of
#'   a line of text. Default is -0.1 for small inward ticks.
#' @param mgp Numeric vector of length 3 for margin line for axis title,
#'   labels and line. Default is c(1.5, 0.35, 0) for compact layout.
#' @param cex Numeric scaling factor for text and symbols. Default is 1.
#' @param cex.lab Numeric scaling factor for axis labels. Default is 9/12
#'   (0.75).
#' @param cex.axis Numeric scaling factor for axis annotation. Default is
#'   7/12 (≈0.58).
#' @param bty Character specifying box type around plot. Default is "n"
#'   (no box).
#' @param lend Integer specifying line end style. Default is 2 (butt line
#'   ends).
#' @param ... Additional parameters passed to par().
#'
#' @return
#' Invisibly returns the previous par settings, allowing for restoration if
#' needed.
#'
#' @details
#' This function implements a consistent and aesthetically pleasing style for
#' R plots with the following key features:
#' * Equal margins on all sides (0.9 inches)
#' * Horizontal axis labels
#' * Small, inward-facing tick marks
#' * Compact axis title and label positioning
#' * Slightly reduced text sizes for better proportions
#' * No box around the plot
#' * Clean line endings
#'
#' The defaults are optimized for 3x3 inch square plots but work well with
#' other dimensions.
#'
#' @examples
#' # Basic usage
#' pdf("example1.pdf", width = 3, height = 3)
#' plot(1:10, main = "Before vl_par()")
#' vl_par()
#' plot(1:10, main = "After vl_par()")
#' dev.off()
#'
#' # Multiple plots with custom margins
#' pdf("example2.pdf", width = 6, height = 6)
#' par(mfrow = c(2,2))
#' vl_par(mai = c(0.5, 0.5, 0.5, 0.5))
#' for(i in 1:4) {
#'   plot(rnorm(100), main = paste("Plot", i))
#' }
#' dev.off()
#'
#' # Save and restore previous parameters
#' old_par <- vl_par()
#' plot(1:10)
#' par(old_par)  # restore original settings
#'
#' @seealso
#' \code{\link{par}} for all available graphical parameters
#'
#' @export
vl_par <- function(mai= c(.9, .9, .9, .9),
                 las= 1,
                 tcl= -0.1,
                 mgp= c(1.5, 0.35, 0),
                 cex= 1,
                 cex.lab= 9/12,
                 cex.axis= 7/12,
                 bty= "n",
                 lend= 2,
                 font.main= 1,
                 ...)
{
  par(mai= mai,
      las= las,
      tcl= tcl,
      mgp= mgp,
      cex= cex,
      cex.lab= cex.lab,
      cex.axis= cex.axis,
      bty= bty,
      lend= lend,
      font.main= font.main,
      ...)
}
