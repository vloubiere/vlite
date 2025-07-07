#' Plots balloons legend
#'
#' @param breaks Numeric values specifying size breaks.
#' @param labels Labels corresponding to each size.
#' @param main Title
#' @param cex cex parameter affecting balloons sizes.
#' @param left left position.
#' @param top top position.
#' @param height height
#' @param ticks.cex cex expansion factor for ticks labels
#' @param main.cex cex expansion factor for legend title
#'
#' @return Plots heatkey
#' @export
balloonskey <- function(breaks,
                        labels,
                        main= NA,
                        cex,
                        left= par("usr")[2],
                        top= par("usr")[4],
                        ticks.cex= 0.6,
                        main.cex= 1)
{
  # Checks
  if(length(breaks) != length(labels))
    stop("breaks and labels should have the same legnth")

  # Apply expansion value to breaks
  breaks <- breaks*cex

  # Compute points width in order to center and space them correctly
  width <- strwidth("M", cex= 0.55)*max(abs(breaks))
  center.top <- top-strheight("M", cex= 0.7)/2*max(abs(breaks), na.rm= TRUE)
  center.height <- strheight("M", cex= 0.7)*max(abs(breaks), na.rm= TRUE)*(length(breaks)-1)

  # Y positions
  y <- seq(center.top-center.height,
           center.top,
           length.out=length(breaks))

  # Plotting points
  points(rep(left+width/2, length(breaks)),
         y,
         pch= ifelse(breaks>=0, 21, 22), # Positive values round, neg are squared
         cex= abs(breaks)+.1, # Absolute pervent neg. values from being removed, avoid 0s
         xpd= T)

  # Add break labels
  text(rep(left+width+strwidth("M", cex= 0.7)/2, length(labels)),
       y,
       labels,
       pos= 4,
       offset= 0.25,
       xpd=T,
       cex= ticks.cex)

  # Title
  text(left,
       top+strheight("M"),
       main,
       pos= 4,
       cex= main.cex,
       offset= 0,
       xpd= T)
}
