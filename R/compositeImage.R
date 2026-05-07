#' Title
#'
#' @param im 
#' @param red 
#' @param green 
#' @param blue 
#' @param white 
#'
#' @returns
#' @export
#'
#' @examples
compositeImage <- function(im, red= NULL, green= NULL, blue= NULL, white= NULL) {
  stopifnot(inherits(im, "AnnotatedImage") | inherits(im, "array"))
  stopifnot(length(dim(im)) == 3)
  
  comp <- EBImage::rgbImage(
    red = if(is.null(red) & is.null(white)) NULL else Reduce(`+`, lapply(unlist(c(red, white)), function(i) im[, , i])),
    green = if(is.null(green) & is.null(white)) NULL else Reduce(`+`, lapply(unlist(c(green, white)), function(i) im[, , i])),
    blue = if(is.null(blue) & is.null(white)) NULL else Reduce(`+`, lapply(unlist(c(blue, white)), function(i) im[, , i]))
  )
  
  return(comp)
}