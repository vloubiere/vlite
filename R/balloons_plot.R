#' Ballons plot
#'
#' Starting from two color.varrices (dot size and colors, respectively), generates a balloon plot
#'
#' @param size.var A color.varrix of numeric values controlling balloons'sizes.
#' @param color.var A color.varrix of numeric values controlling balloons'colors.
#' @param color.breaks A numeric vector specifying the breakpoints for color mapping.
#'   Defaults to 21 evenly spaced breaks spanning the range of `color.var`.
#' @param col A vector of colors corresponding to the values in `breaks`. If `NULL`,
#'   a diverging blue-white-red palette is used for data spanning positive and negative values,
#'   or a sequential viridis-like palette.
#' @param cex Expansion factor for balloons. default= 1.
#' @param main Title. Default= NA
#' @param xaxt See ?par(). Default= "s".
#' @param yaxt See ?par(). Default= "s".
#' @param gap.axis See ?par(). Default= 0
#' @examples
#' color.var <- size.var <- matrix(-3:5, ncol= 3)
#' rownames(size.var) <- c("r1", "r2", "r3")
#' colnames(size.var) <- c("A", "B", "C")
#' balloons_plot(size.var= size.var, color.var= color.var)
#'
#' @return Balloon plot
#' @export
balloons_plot <- function(size.var,
                          color.var,
                          color.breaks= NULL,
                          col= NULL,
                          cex= 1,
                          main= NULL,
                          xlab= NULL,
                          ylab= NULL,
                          gap.axis= 0,
                          legend.left.pos= NULL,
                          size.legend.title= NULL,
                          size.legend.breaks= NULL,
                          color.legend.title= NULL)
{
  # Checks ----
  if(is.null(colnames(size.var))) {
    colnames(size.var) <- if(!is.null(colnames(color.var))) {
      colnames(color.var)
    } else {
      seq(ncol(size.var))
    }
  }

  # Default color.breaks ----
  if(is.null(color.breaks)) {
    color.breaks <- if(min(color.var, na.rm= TRUE)<0 & max(color.var, na.rm= TRUE)>0) {
      lims <- max(abs(color.var), na.rm= TRUE)
      # Centered on 0
      seq(-lims,
          lims,
          length.out= 21)
    } else {
      # Otherwise
      seq(min(color.var, na.rm= TRUE),
          max(color.var, na.rm= TRUE),
          length.out= 21)
    }
  }

  # Default colors ----
  if(is.null(col)) {
    col <- if(color.breaks[1] < 0 & color.breaks[length(color.breaks)] > 0){
      # Positive and negative values
      c("royalblue1", "white", "red")
    } else {
      # Otherwise
      c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF")
    }
  }
  col <- colorRampPalette(col)(length(color.breaks))

  # Init plot ----
  plot(x= col(size.var),
       y= row(size.var),
       type= "n",
       xaxt= "n",
       yaxt= "n",
       xlab= NA,
       ylab= NA,
       frame= FALSE)

  # Lines ----
  segments(1:ncol(size.var),
           1,
           1:ncol(size.var),
           nrow(size.var),
           xpd=T)
  segments(1,
           1:nrow(size.var),
           ncol(size.var),
           1:nrow(size.var),
           xpd=T)

  # Compute balloons positions ----
  x <- unlist(col(size.var))
  y <- unlist(row(size.var)[nrow(size.var):1,])
  # Sizes
  sizes <- unlist(size.var)
  # Compute balloons colors
  Cc <- circlize::colorRamp2(color.breaks, col)
  b.col <- Cc(unlist(color.var))
  b.col[is.na(b.col)] <- "lightgrey"

  # Plot balloons ----
  points(x,
         y,
         col= b.col,
         pch= ifelse(sizes>0, 19, 15),
         cex= abs(sizes)*cex, # Multiply by expansion factor
         xpd= T)

  # Axes ----
  tiltAxis(x = seq(ncol(size.var)),
           labels = colnames(size.var))
  axis(side= 2,
       at= rev(seq(nrow(size.var))),
       labels = rownames(size.var),
       gap.axis= gap.axis,
       lwd= NA)

  # Axes legends ----
  if(!is.null(main))
    title(main= main)
  if(!is.null(xlab))
    title(xlab = xlab)
  if(!is.null(ylab))
    title(ylab = ylab)

  # Caompute default left position for legend ----
  if(is.null(legend.left.pos)) {
    legend.left.pos <- par("usr")[2]+diff(grconvertX(c(0,1), "line", "user"))
  }

  # Add color legend ----
  # breaks have to be centered
  d <- diff(color.breaks)
  centered.breaks <- c(
    color.breaks[1] - d[1]/2,
    (color.breaks[-1] + color.breaks[-length(color.breaks)])/2,
    color.breaks[length(color.breaks)] + d[length(d)]/2
  )
  heatkey(breaks = centered.breaks,
          col = col,
          main = size.legend.title,
          position = "right")

  # Compute size breaks (only affects the legend) ----
  if(is.null(size.legend.breaks)) {
    size.legend.breaks <- range(unlist(size.var), na.rm= TRUE)
    size.legend.breaks <- axisTicks(size.legend.breaks, log= F, nint = 4)
  }

  # Add sizes legend ----
  heatkey.line.width <- 4
  heatkey.top <- par("usr")[4]
  heatkey.bot <- heatkey.top-diff(grconvertY(c(0, heatkey.line.width+2), "line", "user"))
  balloonskey(breaks = size.legend.breaks,
              labels = size.legend.breaks,
              main = color.legend.title,
              cex= cex,
              left= legend.left.pos,
              top= heatkey.bot)
}
