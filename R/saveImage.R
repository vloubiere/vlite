
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
saveImage <- function(
    im,
    jpeg.file,
    scale.metadata= NULL,
    scale.bar= 100,
    horizontal.lines= NULL,
    box= FALSE,
    legend= NULL,
    legend.col= "white",
    cex= 1,
    lwd= 1
) {
  # Apply scaling factor
  cex <- cex*20
  lwd <- lwd*20
  
  # Open JPEG
  par(xaxs= "i", yaxs= "i")
  jpeg(jpeg.file, width = dim(im)[1], height = dim(im)[2])
  
  # Plot image
  plot(im)
  
  # Add scale bar
  if(!is.na(scale.bar)) {
    if(is.null(scale.metadata))
      message("No metadta provided for scale bar. Skipped.") else {
        
        # Import metadata
        meta <- getMetadataImage(scale.metadata)
        
        # Compute padding
        x.adj <- meta$sizeX/20
        y.adj <- meta$sizeY/40
        
        # Add scale bar
        segments(
          x0= meta$sizeX-x.adj,
          y0= meta$sizeY-y.adj,
          x1= meta$sizeX-x.adj-(scale.bar/meta$pxSizeX_um),
          y1= meta$sizeX-y.adj,
          col= "white",
          lwd= lwd,
          lend= 2
        )
      }
  }
  
  # Add lines
  if(!is.null(horizontal.lines)) {
    abline(
      h= horizontal.lines,
      col= "white",
      lwd= lwd
    )
  }
  
  # Box
  if(box) {
    box(lwd= lwd*2, col= "white", which = "figure")
  }
  
  # Add legend
  if(!is.null(legend)) {
    vl_legend(
      x = par("usr")[1]-strwidth("M", cex= cex)*2,
      y = par("usr")[4]+strheight("M", cex= cex),
      legend= legend,
      text.col= legend.col,
      cex= cex
    )
  }
  dev.off()
}
