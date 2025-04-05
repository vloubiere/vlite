#' Ballons plot
#'
#' Starting from two matrices (dot size and colors, respectively), generates a balloon plot
#'
#' @param var1 A matrix of numeric values controling balloons'sizes.
#' @param var2 A matrix of numeric values controling balloons'colors.
#' @param color.breaks Breaks used for coloring (see var 2).
#' @param col Vector of colors used for coloring
#' @param cex Expansion factor for balloons. default= 1.
#' @param main Title. Default= NA
#' @param xaxt See ?par(). Default= "s".
#' @param yaxt See ?par(). Default= "s".
#' @param gap.axis See ?par(). Default= 0
#' @examples
#' var1 <- matrix(-3:5, ncol= 3)
#' var2 <- matrix(-3:5, ncol= 3)
#' balloons_plot(var1= var1,
#'                  var2= var2)
#'
#' @return Balloon plot
#' @export
balloons_plot <- function(var1,
                             var2,
                             color.breaks= NULL,
                             col= c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF"),
                             cex= 1,
                             main= NULL,
                             xlab= NULL,
                             ylab= NULL,
                             gap.axis= 0,
                             var1.legend.title= NULL,
                             var2.legend.title= NULL,
                             legend.left.pos= NULL,
                             legend.size.breaks= NULL)
{
  # Checks
  if(is.null(color.breaks)) {
    color.breaks <- range(var2, na.rm= T)
    if(length(unique(color.breaks))==1)
      color.breaks <- color.breaks+c(-0.5, 0.5)
    color.breaks <- seq(min(color.breaks, na.rm= T),
                        max(color.breaks, na.rm= T),
                        length.out= length(col))
  }
  if(is.null(legend.left.pos)) {
    legend.left.pos <- par("usr")[2]+diff(grconvertX(c(0,1), "line", "user"))
  }
  if(is.null(legend.size.breaks)) {
    legend.size.breaks <- range(unlist(var1), na.rm= TRUE)
    legend.size.breaks <- axisTicks(legend.size.breaks, log= F, nint = 4)
  }
  if(is.null(colnames(var))){
    colnames(var1) <- seq(ncol(var1))
  }

  # Init plot
  plot(x= col(var1),
       y= row(var1),
       type= "n",
       xaxt= "n",
       yaxt= "n",
       xlab= NA,
       ylab= NA,
       frame= FALSE)

  # Lines
  segments(1:ncol(var1),
           1,
           1:ncol(var1),
           nrow(var1),
           xpd=T)
  segments(1,
           1:nrow(var1),
           ncol(var1),
           1:nrow(var1),
           xpd=T)


  # Compute balloons positions
  x <- unlist(col(var1))
  y <- unlist(row(var1)[nrow(var1):1,])
  # Sizes
  sizes <- unlist(var1)
  # Colors
  Cc <- circlize::colorRamp2(color.breaks, col)
  col <- Cc(unlist(var2))
  col[is.na(col)] <- "lightgrey"

  # Plot balloons
  points(x,
         y,
         col= col,
         pch= ifelse(sizes>0, 19, 15),
         cex= abs(sizes)*cex, # Multiply by expansion factor
         xpd= T)

  # Axes
  tiltAxis(x = seq(ncol(var1)),
           labels = colnames(var1))
  axis(side= 2,
       at= seq(nrow(var1)),
       labels = rownames(var1),
       gap.axis= gap.axis,
       lwd= NA)

  # Axes legends
  if(!is.null(main))
    title(main= main)
  if(!is.null(xlab))
    title(xlab = xlab)
  if(!is.null(ylab))
    title(ylab = ylab)

  # Add color legend
  heatkey.line.width <- 4
  heatkey.top <- par("usr")[4]
  heatkey(col = col,
          breaks = color.breaks,
          top = heatkey.top,
          left = legend.left.pos,
          main = var1.legend.title,
          height = heatkey.line.width)

  # Add sizes legend
  heatkey.bot <- heatkey.top-diff(grconvertY(c(0, heatkey.line.width+2), "line", "user"))
  balloonskey(size.breaks = legend.size.breaks,
              cex= cex,
              labels = legend.size.breaks,
              left= legend.left.pos,
              top= heatkey.bot,
              main = var2.legend.title)
}
