#' Title
#'
#' @param im 
#' @param clip 
#'
#' @returns
#' @export
#'
#' @examples
clipImage <- function(im, clip= 0) {
  stopifnot(inherits(im, "AnnotatedImage") | inherits(im, "array"))
  
  # For each channel
  for(i in seq_len(dim(im)[3])) {
    .c <- if(length(dim(im)) == 4)
      im[, , i, ] else
        im[, , i]
    
    q <- quantile(as.numeric(.c), c(clip, 1 - clip))
    .c <- (.c - q[1]) / diff(q)
    .c[.c < 0] <- 0
    .c[.c > 1] <- 1
    
    if(length(dim(im)) == 4)
      im[, , i, ] <- .c else
        im[, , i] <- .c
  }
  
  return(im)
}