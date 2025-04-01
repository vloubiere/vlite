#' Create Tilted Axis Labels
#'
#' @description
#' Creates tilted labels for x-axis tick marks, improving readability when
#' dealing with long or numerous labels. This function provides a more
#' flexible alternative to standard axis labels, particularly useful for
#' heatmaps, barplots, and other visualizations with potentially overlapping
#' labels.
#'
#' @param x Numeric vector specifying x-coordinates for labels.
#' @param y Numeric vector specifying y-coordinates for labels. Defaults to
#'   bottom of plot area adjusted by margin parameters.
#' @param labels Character vector of labels to be plotted.
#' @param srt Numeric angle of rotation in degrees. Default is 45.
#' @param offset Numeric offset from the x position in character widths.
#'   Default is 0.25.
#' @param pos Integer specifying position relative to point:
#'   1=below, 2=left, 3=above, 4=right. Default is 2.
#' @param xpd Logical or NA controlling clipping region. Default is NA
#'   (plotting allowed outside plot region).
#' @param cex Numeric character expansion factor. Default is par("cex.axis").
#' @param ... Additional arguments passed to text().
#'
#' @return
#' Invisibly returns NULL. Called for its side effect of adding tilted
#' labels to a plot.
#'
#' @details
#' The function automatically calculates appropriate y-positions for labels
#' based on the current plot parameters if not specified. It is particularly
#' useful for:
#' * Long category names
#' * Dense x-axis labels
#' * Improving readability in complex plots
#' * Custom axis label positioning
#'
#' @examples
#' # Basic usage with long labels
#' plot(1:5, xaxt = "n", xlab = "")
#' tiltAxis(1:5, labels = paste("Long Label", 1:5))
#'
#' @seealso
#' \code{\link{text}} for underlying text plotting function
#' \code{\link{par}} for graphical parameters
#' \code{\link{gImage}} for heatmap plotting function
#'
#' @export
tiltAxis <- function(x,
                     y= rep(par("usr")[3], length(labels))-diff(grconvertY(c(0, par("mgp")[2]), "line", "user")),
                     labels,
                     srt= 45,
                     offset= 0.25,
                     pos= 2,
                     xpd= NA,
                     cex= par("cex.axis"),
                     ...)
{
  text(x,
       y,
       labels= labels,
       srt= srt,
       offset= -offset,
       pos= pos,
       xpd= xpd,
       cex= cex,
       ...)
}
