#' Creates a composite image by merging different channels
#'
#' @param im The image object from which channels should be extracted.
#' @param red Index of the color chanel that should be shown in red. Default= NULL.
#' @param green Index of the color chanel that should be shown in green.  Default= NULL.
#' @param blue Index of the color chanel that should be shown in blue. Default= NULL.
#' @param white Index of the color chanel that should be shown in white. Default= NULL.
#' If all colors are set to FALSE, then an empty image is returned.
#'
#' @returns
#' @export
#'
#' @examples
compositeImage <- function(im, red= NULL, green= NULL, blue= NULL, white= NULL) {
  # Checks ----
  stopifnot(inherits(im, "AnnotatedImage") | inherits(im, "array"))
  stopifnot(length(dim(im)) == 3)
  missing.channels <- setdiff(unlist(c(red, green, blue, white)), seq(dim(im)[3]))
  if(length(missing.channels))
    stop(paste("Channels", paste0(missing.channels, collapse = ","), "are missing"))
  # Is image empty?
  is.empty <- is.null(red) & is.null(green) & is.null(blue) & is.null(white)
  
  # Assemble image ----
  comp <- if(is.empty) 
    EBImage::rgbImage(
      red = matrix(0, nrow= nrow(im[,,1]), ncol= nrow(im[,,1])),
    ) else 
      EBImage::rgbImage(
        red = if(is.null(red) & is.null(white)) NULL else Reduce(`+`, lapply(unlist(c(red, white)), function(i) im[, , i])),
        green = if(is.null(green) & is.null(white)) NULL else Reduce(`+`, lapply(unlist(c(green, white)), function(i) im[, , i])),
        blue = if(is.null(blue) & is.null(white)) NULL else Reduce(`+`, lapply(unlist(c(blue, white)), function(i) im[, , i]))
      )
  
  return(comp)
}