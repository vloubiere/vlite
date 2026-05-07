#' Title
#'
#' @param file 
#'
#' @returns
#' @export
#'
#' @examples
importImage <- function(file) {
  stopifnot(length(file)==1)
  # Import metadata
  meta <- RBioFormats::read.metadata(file)
  names <- if(grepl(".lif$", file)) {
    make.unique(sapply(meta@.Data, function(x) x$seriesMetadata$`Image name`))
  } else if(grepl(".czi$", file)) {
    gsub("(.*)-ApoTome.*", "\\1", basename(file))
  }
  
  # Import images
  im <- lapply(seq(names), function(i) read.image(file, series = i, normalize = TRUE))
  names(im) <- names
  message(length(im), " image(s) imported into a list.")
  
  # Return
  return(im)
}