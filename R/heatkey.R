#' Add Color Key (Legend) to Heatmap
#'
#' Plots a color key (legend) for heatmaps. This function is designed to work alongside heatmap visualizations
#' and should be called after the main heatmap plot.
#'
#' @param breaks A numeric vector specifying the breakpoints for color mapping.
#' @param col A vector of colors corresponding to the values in `breaks`.
#' @param left Numeric position for the left edge of the color key. Default is calculated based on the current plot region.
#' @param top Numeric position for the top edge of the color key. Default is the top of the current plot region.
#' @param height Numeric height of the color key in lines. Default is `4`.
#' @param width Numeric width of the color key in lines. Default is `0.75`.
#' @param main Character string for the title of the color key. Default is `NA` (no title).
#' @param cex Numeric scaling factor for the size of the color key text. Default is `1`.
#'
#' @return
#' Adds a color key to the current plot. No return value.
#'
#' @details
#' The function creates a vertical color key that maps colors to values in your heatmap. The `breaks` vector defines
#' the midpoints of the ranges for each color in `col`. The color key is customizable in terms of position, size, and appearance.
#'
#' @examples
#' # Create example matrix
#' set.seed(1234)
#' test <- matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] <- test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] <- test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] <- test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) <- paste("Test", 1:10, sep = "")
#' rownames(test) <- paste("Gene", 1:20, sep = "")
#'
#' # Basic heatmap with default color key
#' vl_par()
#' vl_image(test)
#' heatkey(col = heat.colors(100), breaks = seq(-3, 3, length.out = 100))
#'
#' @seealso
#' \code{\link{vl_image}} for creating the main heatmap.
#'
#' @export
heatkey <- function(breaks,
                    col,
                    left= par("usr")[2]+diff(grconvertX(c(0,1), "line", "user")),
                    top= par("usr")[4],
                    width= 0.75,
                    height= 4,
                    cex= 1,
                    main= NA)
{
  # Compute plotting positions ----
  line.height <- diff(grconvertY(c(0, height), "line", "user"))
  ybottom <- top-line.height
  ypos <- seq(ybottom, top, length.out= length(breaks)+1)
  line.width <- diff(grconvertX(c(0, width), "line", "user"))
  xright <- left+line.width

  # Plot key ----
  rect(xleft = left,
       ybottom = ypos[-length(ypos)],
       xright = xright,
       ytop = ypos[-1],
       xpd= NA,
       col= col,
       border= NA)

  # Add border ----
  rect(xleft = left,
       ybottom = ybottom,
       xright = xright,
       ytop = top,
       xpd= NA)

  # Compute ticks ----
  ticks <- axisTicks(range(breaks),
                     log= F,
                     nint = 3)

  # Compute center y position of highest and lowest breaks rectangles ----
  min.pos <- mean(data.table::first(ypos, n = 2)) # Lowest break
  max.pos <- mean(data.table::last(ypos, n = 2)) # Highest break
  span <- max.pos-min.pos # Span bewteen lowest and highest
  # Compute ticks y posistions
  ypos.t <- min.pos+span*((ticks-breaks[1])/diff(range(breaks)))

  # Add labels ----
  text(x= left+line.width,
       y= ypos.t,
       labels = ticks,
       cex= cex*0.7,
       pos= 4,
       xpd= NA,
       offset= 0.25)
  segments(x0 = left+line.width,
           y0 = ypos.t,
           x1 = left+line.width+line.width/8,
           y1 = ypos.t,
           xpd= T)
  text(x = left,
       y = top+strheight(main),
       labels = main,
       pos = 4,
       cex = cex,
       xpd = NA,
       offset = 0)
}
