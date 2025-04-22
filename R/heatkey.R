#' Add Color Key (Legend) to Heatmap
#'
#' Plots a color key (legend) for heatmaps. This function is designed to work alongside heatmap visualizations
#' and should be called after the main heatmap plot.
#'
#' @param breaks A numeric vector specifying the breakpoints for color mapping.
#' @param col A vector of colors corresponding to the values in `breaks`.
#' @param pos Position of the heatkey. Can either be "right" or "top". Default= "right".
#' @param adj.x If specified, x plotting positions are adjusted by adj.x * line.width Default= 0.
#' @param adj.y If specified, y plotting positions are adjusted by adj.y * line height. Default= 0.
#' @param thickness Numeric width of the color key in lines. Default= 0.75.
#' @param length Numeric height of the color key in lines. Default= 4.
#' @param main Character string for the title of the color key. Default= NA.
#' @param compute.ticks Should nice ticks be computed? By default, they will if provided breaks are equally spaced.
#' @param cex Numeric scaling factor for the size of the color key text. Default= 1.
#'
#' @return
#' Adds a color key to the current plot, at one line of distance from the plot border. No return value.
#'
#' @details
#' The function creates a vertical color key that maps colors to values in your heatmap. The breaks vector defines
#' the midpoints of the ranges for each color in col. The color key is customizable in terms of position, size, and appearance.
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
                    position= "right",
                    adj.x= 0,
                    adj.y= 0,
                    thickness= 0.75,
                    length= 4,
                    main= NA,
                    compute.ticks,
                    cex= 1)
{
  # Checks ----
  if(!position %in% c("right", "top"))
    stop("position should be one of 'right' or 'top'.")

  # Compute normalized breaks ----
  min.break <- min(diff(breaks))
  span <- diff(range(breaks))
  norm.breaks <- seq(min(breaks),
                     max(breaks),
                     length.out= round(span/min.break)+1)

  # Compute normalized colors ----
  col <- cut(norm.breaks,
             c(-1, breaks),
             include.lowest = T,
             labels = col)
  col <- as.character(col)

  # Compute line width and height (used as reference) ----
  line.width <- diff(grconvertX(c(0, 1), "line", "user"))*cex
  line.height <- diff(grconvertY(c(0, 1), "line", "user"))*cex

  # Compute plotting positions ----
  if(position=="top") {
    # X pos
    xleft <- mean(par("usr")[1:2])-(length/2*line.width)
    xright <- xleft+length*line.width
    pos <- seq(xleft,
               xright,
               length.out= length(norm.breaks)+1)
    xleft <- pos[-length(pos)]
    xright <- pos[-1]
    # Y pos
    ybottom <- par("usr")[4]+line.height
    ytop <- ybottom+thickness*line.height
  } else if(position=="right") {
    # X pos
    xleft <- par("usr")[2]+line.width
    xright <- xleft+thickness*line.width
    # Y pos
    ytop <- par("usr")[4]
    ybottom <- ytop-length*line.height
    pos <- seq(ybottom,
               ytop,
               length.out= length(norm.breaks)+1)
    ybottom <- pos[-length(pos)]
    ytop <- pos[-1]
  }

  # Adjust position ----
  xadj <- adj.x*line.width
  yadj <- adj.y*line.height

  # Plot key ----
  rect(xleft = xleft+xadj,
       ybottom = ybottom+yadj,
       xright = xright+xadj,
       ytop = ytop+yadj,
       xpd= NA,
       col= col,
       border= NA)

  # Add border ----
  rect(xleft = xleft[1]+xadj,
       ybottom = ybottom[1]+yadj,
       xright = rev(xright)[1]+xadj,
       ytop = rev(ytop)[1]+yadj,
       xpd= NA)

  # Add title ----
  text(x = ifelse(position=="top", mean(pos), xleft[1])+xadj,
       y = rev(ytop)[1]+ifelse(position=="top", line.height, line.height/2)+yadj,
       labels = main,
       pos = ifelse(position=="top", 3, 4),
       cex = cex,
       xpd = NA,
       offset = 0)

  # Check if nice ticks should be computed ----
  if(missing(compute.ticks)) {
    tol <- diff(range(breaks))/1000
    diffs <- diff(breaks)
    equally.spaced <- all.equal(target = diffs,
                                current = rep(diffs[1], length(diffs)),
                                tolerance = tol)
    compute.ticks <- isTRUE(equally.spaced)
  }

  # Compute ticks ----
  ticks <- if(compute.ticks) {
    axisTicks(range(breaks),
              log= F,
              nint = 3)
  } else {
    breaks
  }

  # Compute tick positions ----
  min.pos <- mean(pos[1:2]) # Lowest break
  max.pos <- mean(rev(pos)[1:2]) # Highest break
  pos.t <- ticks/span*(max.pos-min.pos)+min.pos

  # Plot ticks ----
  if(position=="top")
  {
    x0 <- x1 <- pos.t
    y0 <- ytop
    y1 <- ytop+line.height/8
  } else if(position=="right") {
    x0 <- xright
    x1 <- xright+line.width/8
    y0 <- y1 <- pos.t
  }

  # Add ticks ----
  text(x= x1+xadj,
       y= y0+yadj,
       labels = ticks,
       cex= cex*0.7,
       pos= ifelse(position=="top", 3, 4),
       xpd= NA,
       offset= 0.25)
  segments(x0 = x0+xadj,
           y0 = y0+yadj,
           x1 = x1+xadj,
           y1 = y1+yadj,
           xpd= T)
}
