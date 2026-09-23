#' Title
#'
#' @param im Two-dimensional image to plot.
#' @param pxSizeX_um Pixel width in micrometres.
#' @param scale.bar Scale-bar length in micrometres. Use `NA` or `NULL` to omit it.
#'
#' @returns
#' @export
#'
#' @examples
plotImage <- function(
    im,
    pxSizeX_um= NULL,
    scale.bar= 100
) {
  stopifnot(inherits(im, c("matrix", "array", "Image")))
  
  # Adjust plotting parameters
  if(nrow(im)==ncol(im))
    par(xaxs= "i", yaxs= "i", pty= "s") else
      par(xaxs= "i", yaxs= "i")
  
  # Plot image
  EBImage::display(
    im,
    method= "raster",
    all= all # all frames of a stacked EBI image will be arranged in a grid
  )
  
  # Add scale bar
  if(!is.null(scale.bar) && !is.na(scale.bar)) {
    if(is.null(pxSizeX_um)) {
      warning("pxSizeX_um missing -> scale bar skipped!")
    } else {
      right <- par("usr")[2]-diff(par("usr")[c(1,2)])/20
      y <- par("usr")[3]+diff(par("usr")[c(3,4)])/40
      segments(
        x0= right,
        y0= y,
        x1= right-(scale.bar/pxSizeX_um),
        y1= y,
        col= "white",
        lwd= 2,
        lend= 2
      )
    }
  }
}
