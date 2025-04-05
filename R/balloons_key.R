#' Plots balloons legend
#'
#' @param size.breaks Balloons sizes.
#' @param cex cex parameter affecting balloons sizes.
#' @param labels Labels corresponding to each size.
#' @param left left position.
#' @param top top position.
#' @param height height
#' @param main Title
#' @param ticks.cex cex expansion factor for ticks labels
#' @param main.cex cex expansion factor for legend title
#'
#' @return Plots heatkey
#' @export
balloonskey <- function(size.breaks,
                        cex,
                        labels,
                        left= par("usr")[2],
                        top= par("usr")[4],
                        ticks.cex= 0.6,
                        main.cex= 1,
                        main= NA)
{
  # Checks
  if(length(size.breaks) != length(labels))
    stop("size.breaks and labels should have the same legnth")

  # Apply expansion value to size.breaks
  size.breaks <- size.breaks*cex

  # Compute points width in order to center and space them correctly
  width <- strwidth("M", cex= 0.55)*max(abs(size.breaks))
  center.top <- top-strheight("M", cex= 0.7)/2*max(abs(size.breaks), na.rm= TRUE)
  center.height <- strheight("M", cex= 0.7)*max(abs(size.breaks), na.rm= TRUE)*(length(size.breaks)-1)

  # Y positions
  y <- seq(center.top-center.height,
           center.top,
           length.out=length(size.breaks))

  # Plotting points
  points(rep(left+width/2, length(size.breaks)),
         y,
         cex= abs(size.breaks)+0.1,
         pch= ifelse(size.breaks>=0, 21, 22),
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
