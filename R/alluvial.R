#' Create an Alluvial (Sankey) Diagram
#'
#' @description
#' Creates an alluvial (Sankey) diagram showing the flow between categories across multiple variables.
#' The diagram consists of vertical bars representing categories and connecting polygons showing the
#' relationships between consecutive categories.
#'
#' @param x A data.table where each column represents a categorical variable and each row is an entry (see examples).
#' @param bars.widths Numeric vector specifying the width of the vertical bars. If a single value is provided,
#'        it will be recycled for all columns. Default is 0.3
#' @param col Vector of colors used for the color gradient of categories. Default is rainbow(12)
#' @param xlab Character string for x-axis label. Default is "Categories"
#' @param ylab Character string for y-axis label. Default is "N"
#' @param keep.order Logical indicating whether to maintain the original order of factor levels. Default is FALSE
#' @param border Color of the border for vertical bars. Default is "black"
#' @param xlim Numeric vector of length 2 giving the x coordinate range. Default is NULL (automatically computed)
#' @param ylim Numeric vector of length 2 giving the y coordinate range. Default is NULL (automatically computed)
#' @param xaxt Character specifying x-axis type. Default is "n" (no axis)
#' @param yaxt Character specifying y-axis type. Default is "s" (standard axis)
#' @param show.labels Logical indicating whether to show category labels. Default is TRUE
#' @param labels.min.N Numeric threshold for showing labels (only shows labels for categories with N >= labels.min.N).
#'        Default is 0
#' @param labels.cex Numeric scaling factor for label size. Default is 0.8
#' @param alpha.categories Numeric value between 0 and 1 for the transparency of category bars. Default is 0.8
#' @param alpha.connections Numeric value between 0 and 1 for the transparency of connecting polygons. Default is 0.5
#' @param space Numeric value specifying the vertical space between category groups. Default is 0
#'
#' @return
#' Creates a plot on the current graphics device. No value is returned.
#'
#' @details
#' The function creates an alluvial diagram where vertical bars represent categories and connecting
#' polygons show the flow between categories. The width of the bars and connections is proportional
#' to the number of observations in each category or transition.
#'
#' The diagram can be customized using various parameters for colors, spacing, labels, and transparency.
#' The layout is automatically computed to optimize the use of space and visibility of connections.
#'
#' @examples
#' test <- data.table(c("A", "A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
#'                    c("B", "A", "A", "B", "C", "B", "B", "C", "C", "A", "A"),
#'                    c("A", "A", "B", "B", "B", "B", "B", "B", "B", "C", "A"))
#' alluvial(test)
#' alluvial(test, bars.widths = .1, space = 1)
#'
#' @import data.table
#' @importFrom graphics plot rect polygon text axis
#' @importFrom grDevices adjustcolor colorRampPalette rainbow
#'
#' @seealso
#' \code{\link[graphics]{plot}}, \code{\link[graphics]{polygon}}
#'
#' @note
#' The input data.table is copied internally to prevent modification of the original data.
#' The function works best with categorical data with a reasonable number of categories.
#'
#' @author Your Name
#'
#' @export
alluvial <- function(x,
                     bars.widths= .3,
                     col= rainbow(12),
                     xlab= "Categories",
                     ylab= "N",
                     keep.order= FALSE,
                     border= "black",
                     xlim= NULL,
                     ylim= NULL,
                     xaxt= "s",
                     yaxt= "s",
                     show.labels= TRUE,
                     labels.pos= "center",
                     labels.min.N= 0,
                     labels.cex= .8,
                     alpha.categories= .8,
                     alpha.connections= .5,
                     space= 0)
{
  # Copy for encapsulation
  dat <- data.table::copy(x)

  # Check ordering
  if(keep.order)
    dat[, (names(dat)) := lapply(.SD, function(x) factor(x, unique(x)))]
  setorderv(dat, names(dat))

  # Checks
  if(length(bars.widths)!=ncol(dat))
    bars.widths <- rep(bars.widths, length.out= ncol(dat))
  if(is.null(xlim))
    xlim <- c(1-data.table::first(bars.widths),
              ncol(dat)+data.table::last(bars.widths))
  if(is.null(ylim))
    ylim <- c(0,
              nrow(dat)+max(sapply(dat, data.table::uniqueN)-1)*space)

  # Initiate plot
  plot(NA,
       xlim= xlim,
       ylim= ylim,
       frame= FALSE,
       yaxt= yaxt,
       xaxt= "n",
       xlab= xlab,
       ylab= ylab)
  if(xaxt!="n") {
    axis(1,
         at = seq(dat),
         labels = names(dat),
         lwd = 0,
         padj= -1.25)
  }

  # Initiate labels (plotted last)
  labs <- data.table(x= numeric(),
                     y= numeric(),
                     N= numeric(),
                     lab= character())

  # For each column
  for(i in 1:ncol(dat)) {
    # Make 2-columns object
    .c <- if(i<ncol(dat)) {
      dat[, c(i, i+1), with= FALSE]
    } else
      dat[, c(i, i), with= FALSE]
    setnames(.c, c("V1", "V2"))
    setorderv(.c, c("V1", "V2"), c(-1, -1))

    # Compute colors based on 1st category
    .c[, col:= colorRampPalette(col)(.NGRP)[.GRP], V1]

    # Compute boxes limits coordinates
    .c[, top0:= max(.I)+(.GRP-1)*space, V1]
    .c[, top0:= top0+(ylim[2]-max(top0))/2]
    .c[, bot0:= top0-.N, V1]

    # Compute left connections top coordinates
    .c[, top1:= max(.I), .(V1, V2)]
    .c[, top1:= top1+(.GRP-1)*space, V1]
    .c[, top1:= top1+(ylim[2]-max(top1))/2]

    # Compute right connections top coordinates
    setorderv(.c, "V2", -1)
    .c[, top2:= max(.I), .(V2, V1)]
    .c[, top2:= top2+(.GRP-1)*space, V2]
    .c[, top2:= top2+(ylim[2]-max(top2))/2]

    # Compute connection height
    .c[, h:= .N, .(V1, V2)]

    # Plot Categories
    .c[, {
      rect(xleft= i-bars.widths[i],
           ybottom= bot0[1],
           xright= i+bars.widths[i],
           ytop= top0[1],
           border= border,
           col= adjustcolor(col[1], alpha.f = alpha.categories))
    }, .(top0, bot0, col)]

    # Retrieve labels
    .l <- .c[, {
      .(x= ifelse(labels.pos=="right", i+bars.widths[i], i),
        y= mean(c(top0, bot0)),
        N= .N)
    }, .(lab= V1, top0, bot0)]
    .l$top0 <- .l$bot0 <- NULL
    labs <- rbind(labs, .l)

    # Plot connections
    if(i<ncol(dat)) {
      .c[, {
        # Define control points for bezier curves
        x <- c(i+bars.widths[i], i+1-bars.widths[i+1])
        x <- c(x[1], x[1]+diff(x)/3, x[1]+2*diff(x)/3, x[2])
        y <- c(top1[1], top1[1], top2[1], top2[1])
        cp <- matrix(c(x,y), ncol= 2)
        curve <- bezier::bezier(t = seq(0, 1, length.out = 20), p = cp)
        # Plot connections polygon
        polygon(c(curve[,1], rev(curve[,1])),
                c(curve[,2], rev(curve[,2])-h),
                border= NA,
                col= adjustcolor(col[1], alpha.f = alpha.connections))
      }, .(top1, top2, h, col)]
    }
  }

  # Plot labels
  if(show.labels && any(labs$N>=labels.min.N)) {
    labs[N>=labels.min.N, {
      if(labels.pos=="right") {
        text(x,
             y,
             lab,
             cex= labels.cex,
             pos = 4,
             xpd= TRUE)
      } else {
        text(x,
             y,
             lab,
             cex= labels.cex)
      }
    }]
  }
}
