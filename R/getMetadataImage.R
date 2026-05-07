#' Title
#'
#' @param file 
#'
#' @returns
#' @export
#'
#' @examples
getMetadataImage <- function(file) {
  stopifnot(length(file)==1)
  
  # Import metadata
  meta <- RBioFormats::read.metadata(file)
  
  # Format
  meta <- cbind(as.data.table(meta$coreMetadata), as.data.table(meta$globalMetadata))
  
  # Compute pixel size
  pixelX <- as.numeric(strsplit(meta$`ImageScaling|ImagePixelSize`, ",")[[1]][1])
  pixelY <- as.numeric(strsplit(meta$`ImageScaling|ImagePixelSize`, ",")[[1]][2])
  obj <- as.numeric(meta$`Information|Instrument|Objective|NominalMagnification`)
  magnif <- as.numeric(meta$`Scaling|AutoScaling|CameraAdapterMagnification`)
  magnif <- (magnif*obj)
  pxSizeX <- pixelX/magnif
  pxSizeY <- pixelY/magnif
  
  # Select information
  meta <- meta[, .(sizeX, sizeY, sizeZ, sizeC, magnification= magnif, pxSizeX_um= pxSizeX, pxSizeY_um= pxSizeY)]
  
  # Return
  return(meta)
}
