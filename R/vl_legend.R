#' Add Enhanced Legend to Plot
#'
#' @description
#' A wrapper around base R's legend() function with improved default settings
#' for scientific visualization. Automatically positions the legend in the
#' top-right corner of the plot with clean, minimal styling.
#'
#' @param x Numeric x-coordinate for legend position. Default is right edge
#'   of plot (par("usr")[2]).
#' @param y Numeric y-coordinate for legend position. Default is top edge
#'   of plot (par("usr")[4]).
#' @param legend Character vector of legend labels.
#' @param fill Vector of colors for legend key boxes.
#' @param bty Character specifying box type around legend. Default is "n"
#'   (no box).
#' @param xpd Logical controlling clipping. Default is TRUE (allow legend
#'   outside plot region).
#' @param cex Numeric scaling factor for legend text. Default is 7/12
#'   (≈0.58).
#' @param border Color for borders of legend boxes. Default is NA (no
#'   borders).
#' @param ... Additional arguments passed to legend().
#'
#' @return
#' Invisibly returns a list with legend positioning information (same as
#' legend()).
#'
#' @details
#' This function implements a consistent and clean legend style with:
#' * No box around the legend
#' * No borders around color boxes
#' * Slightly reduced text size
#' * Automatic positioning at plot edge
#' * Ability to extend outside plot region
#'
#' @examples
#' # Basic plot with legend
#' plot(1:10, col = rainbow(10), pch = 16)
#' vl_legend(legend = paste("Group", 1:10),
#'         fill = rainbow(10))
#'
#' @seealso
#' \code{\link{legend}} for the underlying legend function
#' \code{\link{par}} for graphical parameters
#'
#' @export
vl_legend <- function(x= par("usr")[2],
                    y= par("usr")[4],
                    legend,
                    fill,
                    bty= "n",
                    xpd= T,
                    cex= 7/12,
                    border= NA,
                    ...)
{
  legend(x= x,
         y= y,
         legend= legend,
         fill= fill,
         bty= bty,
         xpd= xpd,
         cex= cex,
         border= border,
         ...)
}
