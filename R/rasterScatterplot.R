#' Title
#'
#' @param x
#' @param y default= NULL
#' @param type default= "p"
#' @param frame Default= F
#' @param size Size of the png image. Default= 2000
#' @param res Nominal resolution of png image. default= 600
#' @param ... Extra plotting parameters
#'
#' @return Raster scatterplot
#' @export
#'
#' @examples
#' rasterScatterplot(1:3)
rasterScatterplot <- function(x,
                              y= NULL,
                              pch= 16,
                              col= adjustcolor("grey", .7),
                              cex= .8,
                              type= "p",
                              frame= F,
                              res= 600L,
                              size= 2000L,
                              xlab= NULL,
                              ylab= NULL,
                              xaxt= "s",
                              yaxt= "s",
                              xaxs= "r",
                              yaxs= "r",
                              add= F,
                              ...)
{
  # Initialize plot ----
  if(!add) {
    plot(x= x,
         y= y,
         frame= frame,
         type= "n",
         xlab= xlab,
         ylab= ylab,
         xaxt= "n",
         yaxt= yaxt,
         xaxs= xaxs,
         yaxs= yaxs,
         ...)
    if(xaxt!="n")
      axis(1, padj= -1.25)
  }


  # Extract plot area in both user and physical coordinates ----
  coords <- par("usr")
  gx <- grconvertX(coords[1:2], "user", "inches")
  gy <- grconvertY(coords[3:4], "user", "inches")
  width <- diff(gx)
  height <- diff(gy)
  ratio <- round(c(width, height)/max(c(width, height))*size)

  # Save as png ----
  tmp <- tempfile(fileext = "png")
  png(tmp,
      width = ratio[1],
      height = ratio[2],
      units = "px",
      res = res,
      type="cairo",
      bg = "transparent")
  par(mar = c(0,0,0,0))
  plot.new()
  plot(x= x,
       y= y,
       xlim= coords[1:2],
       ylim= coords[3:4],
       pch= pch,
       col= col,
       cex= cex,
       type= type,
       axes= FALSE,
       xaxs = "i",
       yaxs = "i",
       frame= F)
  dev.off()

  # Plot png ----
  panel <- png::readPNG(tmp)
  rasterImage(panel,
              coords[1],
              coords[3],
              coords[2],
              coords[4])
}
