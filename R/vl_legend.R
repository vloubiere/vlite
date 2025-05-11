#' Add Enhanced Legend to Plot
#'
#' @description
#' A wrapper around base R's legend() function with improved default settings
#' for scientific visualization. Automatically positions the legend in the
#' top-right corner of the plot with clean, minimal styling.
#'
#' @param x Numeric x-coordinate for legend position. Default= par("usr")[2].
#' @param y Numeric y-coordinate for legend position. Default= par("usr")[4].
#' @param legend Character vector of legend labels.
#' @param fill Vector of colors for legend key boxes.
#' @param x.adj If specified, x plotting position is shifted by x.adj*line.width. Default= 0.
#' @param y.adj If specified, y plotting positions is shifted by y.adj*line.height. Default= 0.
#' @param x.log Should the x value be logged? Default= par("xlog").
#' @param y.log Should the x value be logged? Default= par("ylog").
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
vl_legend <- function(x,
                      y,
                      legend,
                      fill,
                      x.adj= 0,
                      y.adj= 0,
                      x.log= par("xlog"),
                      y.log= par("ylog"),
                      bty= "n",
                      xpd= T,
                      cex= 7/12,
                      border= NA,
                      ...)
{
  # Default position ----
  if(missing(x))
    x <- par("usr")[2]
  if(is.numeric(x) && missing(y))
    y <- par("usr")[4]

  # Adjust x and y positions ----
  if(is.numeric(x)) {
    x <- x+diff(grconvertX(c(0,1), "line", "user"))*x.adj
    y <- y+diff(grconvertY(c(0,1), "line", "user"))*y.adj
  }

  # Log if necessary ----
  if(is.numeric(x)) {
    if(x.log)
      x <- 10^x
    if(y.log)
      y <- 10^y
  }

  # Plot legend ----
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
