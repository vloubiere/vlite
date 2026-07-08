#' Title
#'
#' @param path 
#'
#' @returns
#' @export
#'
#' @examples
importImage <- function(path, serie= 1) {
  stopifnot(length(serie)==1)
  
  # Import image
  im <- RBioFormats::read.image(path, series = serie, normalize = TRUE, read.metadata = FALSE)
  im <- im@.Data
  
  # Return
  return(im)
}