#' Title
#'
#' @param im 
#' @param jpeg.file 
#' @param main 
#' @param col 
#' @param cex 
#'
#' @returns
#' @export
#'
#' @examples
plotImage <- function(
    im,
    pxSizeX_um= NULL,
    scale.bar= 100,
    all= FALSE
) {
  stopifnot(class(im)[1] %in% c("matrix", "array", "Image"))
  
  # Adjust plotting parameters
  par(xaxs= "i", yaxs= "i", pty= "s")
  
  # Plot image
  EBImage::display(im, method= "raster", all= all)
  
  # Add scale bar
  if(!is.na(scale.bar)) {
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
