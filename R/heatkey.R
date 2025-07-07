#' Add Color Key (Legend) to Heatmap
#'
#' Plots a color key (legend) for heatmaps. This function is designed to work alongside heatmap visualizations
#' and should be called after the main heatmap plot.
#'
#' @param breaks A numeric vector specifying the breakpoints for color mapping.
#' @param col A vector of colors, which can either be the length of breaks (in which case breaks will be centered)
#'    or have one value less than breaks (in which case breaks will be used as edges).
#' @param labels Labels to be added to the color key. By default, elegant numeric ticks spanning the values in breaks are computed.
#' If numeric labels are provided, they are plotted at the corresponding positions.
#' If a character vector (or factor levels) of length equal to col is provided, labels are centered between each breaks' intervals.
#' @param position Position of the color key. Either "right" (vertical legend) or "top" (horizontal legend). Default= "right".
#' @param adj.x If specified, x plotting positions are adjusted by adj.x*line.width. Default= 0.
#' @param adj.y If specified, y plotting positions are adjusted by adj.y*line.height. Default= 0.
#' @param thickness Numeric width of the color key in lines. Default= 0.75.
#' @param length Numeric height (if vertical) or width (if horizontal) of the color key in lines. Default= 4.
#' @param main Character string for the title of the color key. Default= NA.
#' @param cex Numeric scaling factor for the size of the color key text. Default= 1.
#'
#' @return
#' Adds a color key to the current plot, at one line of distance from the plot border. No return value.
#'
#' @details
#' The function creates a color key (legend) that maps colors to value intervals in your heatmap.
#' The number of colors must be exactly one less than the number of breaks.
#' Labels can be numeric (placed at the corresponding value) or character (placed at the center of each color interval).
#'
#' @examples
#' # Example with the image function
#' image(matrix(1:3, ncol= 1), col= c("blue", "white", "red"))
#'
#' # By default, in this function, breaks are used as edges:
#' heatkey(0:3, c("blue", "white", "red"), main = "test")
#'
#' # But in other situations (for example to emphasize 0 as the middle value), one might prefer to use centered breaks:
#' breaks <- c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5)
#' col <- c("blue", "cornflowerblue", "white", "tomato", "red")
#' image(matrix(-2:2, ncol= 1), breaks= breaks, col= col)
#' heatkey(breaks= breaks, col= col, pos= "top", main= "Values")
#'
#' # Character also be used for character labels
#' heatkey(breaks= 1:3, labels= c("A", "B", "C"), col= c("green", "yellow", "purple"), main= "Test")
#'
#' @export
heatkey <- function(breaks,
                    col,
                    labels= NULL,
                    position= "right",
                    adj.x= 0,
                    adj.y= 0,
                    thickness= 0.75,
                    length= 4,
                    main= NA,
                    cex= 1)
{
  # Checks ----
  if(!position %in% c("right", "top"))
    stop("position should be one of 'right' or 'top'.")
  if(!is.numeric(breaks))
    stop("Breaks should be numeric")
  if(length(breaks)==length(col)) { # Center breaks if necessary
    d <- diff(breaks)
    breaks <- c(
      breaks[1] - d[1]/2,
      (breaks[-1] + breaks[-length(breaks)])/2,
      breaks[length(breaks)] + d[length(d)]/2
    )
  }
  if(length(breaks)!=length(col)+1)
    stop("heatkey: breaks should be the same length as col or one value longer.")
  if(is.null(labels))
    labels <- axisTicks(range(breaks),
                        log= F,
                        nint = 4)
  if(length(labels)!=length(col) && !is.numeric(labels))
    stop("labels should either be the length of col vector or be numeric.")
  if(is.numeric(labels))
    labels <- labels[labels>=min(breaks) & labels<=max(breaks)]

  # Compute line width and heigh (used as reference) ----
  line.width <- diff(grconvertX(c(0, 1), "line", "user"))
  line.height <- diff(grconvertY(c(0, 1), "line", "user"))

  # Compute plotting positions ----
  if(position=="top") {
    # X pos
    xleft <- mean(par("usr")[1:2])-(length/2*line.width*cex)
    xright <- xleft+length*line.width*cex
    pos <- (breaks-min(breaks))/diff(range(breaks))*(xright-xleft)+xleft
    xleft <- pos[-length(pos)]
    xright <- pos[-1]
    # Y pos
    ybottom <- par("usr")[4]+line.height
    ytop <- ybottom+thickness*line.height*cex
  } else if(position=="right") {
    # X pos
    xleft <- par("usr")[2]+line.width
    xright <- xleft+thickness*line.width*cex
    # Y pos
    ytop <- par("usr")[4]
    ybottom <- ytop-length*line.height*cex
    pos <- (breaks-min(breaks))/diff(range(breaks))*(ytop-ybottom)+ybottom
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
       y = rev(ytop)[1]+ifelse(position=="top",
                               line.height*0.75*cex,
                               line.height*cex/2)+yadj,
       labels = main,
       pos = ifelse(position=="top", 3, 4),
       cex = cex,
       xpd = NA,
       offset = 0)

  # Compute ticks ----
  span <- pos[length(pos)]-pos[1]
  pos.t <- if(is.numeric(labels)) {
    (labels-min(breaks))/diff(range(breaks))*span+pos[1]
  } else if(length(labels)==length(col)) {
    pos[-1]-diff(pos)/2
  }

  # Plot ticks ----
  if(position=="top")
  {
    x0 <- x1 <- pos.t
    y0 <- ytop
    y1 <- ytop+line.height*cex/8
  } else if(position=="right") {
    x0 <- xright
    x1 <- xright+line.width*cex/8
    y0 <- y1 <- pos.t
  }

  # Add ticks ----
  text(x= x1+xadj,
       y= y0+yadj,
       labels = labels,
       cex= cex*0.7,
       pos= ifelse(position=="top", 3, 4),
       xpd= NA,
       offset= 0.25*cex)
  segments(x0 = x0+xadj,
           y0 = y0+yadj,
           x1 = x1+xadj,
           y1 = y1+yadj,
           xpd= NA)
}
