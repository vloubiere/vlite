#' Add Color Key (Legend) to Heatmap
#'
#' @description
#' Plots a color key (legend) for heatmaps.
#' This function is designed to work alongside heatmapvisualizations
#' and should be called after the main heatmap plot.
#'
#' @param col Color vector of length(breaks)-1.
#' @param breaks Numeric vector of break points for color mapping, of length(col)+1.
#' @param left Numeric position for the left edge of the color key. Default is
#'   calculated based on the current plot region plus one line width.
#' @param top Numeric position for the top edge of the color key. Default is
#'   the top of the current plot region.
#' @param height Numeric height of the color key in lines (4 by default).
#' @param width Numeric width of the color key in lines (0.75 by default).
#' @param main Character string for the title of the color key (NA by default).
#' @param cex Numeric scaling factor for the heatkey. Default= 1.
#' @param border Color for the border of the color key. Use NA for no border
#'   ("black" by default).
#' @param show.breaks Logical; whether to show break values (TRUE by default).
#'
#' @return
#' Adds a color key to the current plot. No return value.
#'
#' @details
#' The function creates a vertical color key that maps colors to values in your
#' heatmap. It supports both continuous gradients and discrete color categories,
#' and can be customized in terms of position, size, and appearance.
#'
#' @examples
#' # Create example matrix
#' set.seed(1234)
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#'
#' # Basic heatmap with default color key
#' gPar()
#' gImage(test)
#' heatkey()
#'
#' @seealso
#' \code{\link{gHeatmap}} for creating the main heatmap
#'
#' @export
heatkey <- function(col,
                    breaks,
                    left= par("usr")[2]+diff(grconvertX(c(0,1), "line", "user")),
                    top= par("usr")[4],
                    width= 0.75,
                    height= 4,
                    cex= 1,
                    main= NA)
{
  # Compute line width, height adn plotting positions ----
  line.width <- diff(grconvertX(c(0, width), "line", "user"))
  line.height <- diff(grconvertY(c(0, height), "line", "user"))
  ybottom <- top-line.height
  xright <- left+line.width

  # Plot key ----
  rasterImage(image = matrix(rev(col), ncol= 1),
              xleft = left,
              ybottom = ybottom,
              xright = xright,
              ytop = top,
              interpolate = FALSE,
              xpd= NA)

  # Add border ----
  rect(xleft = left,
       ybottom = ybottom,
       xright = xright,
       ytop = top,
       xpd= NA)

  # Compute ticks ----
  ticks <- axisTicks(range(breaks),
                     log= F,
                     nint = 4)

  # Compute y ticks position ----
  height <- diff(c(ybottom, top))
  ypos <- ybottom+height*((ticks-min(breaks))/diff(range(breaks)))

  # Add labels ----
  text(x= left+line.width,
       y= ypos,
       labels = ticks,
       cex= cex*0.7,
       pos= 4,
       xpd= NA,
       offset= 0.25)
  segments(x0 = left+line.width,
           y0 = ypos,
           x1 = left+line.width+line.width/8,
           y1 = ypos,
           xpd= T)
  text(x = left,
       y = top+strheight(main),
       labels = main,
       pos = 4,
       cex = cex,
       xpd = NA,
       offset = 0)
}
