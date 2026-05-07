#' Title
#'
#' @param im 
#' @param start 
#' @param end 
#' @param FUN 
#'
#' @returns
#' @export
#'
#' @examples
maxProjectionImage <- function(im, start = 1, end = dim(im)[4], FUN = pmax) {
  stopifnot(inherits(im, "AnnotatedImage") | inherits(im, "array"))
  
  # If several z-sections
  if(length(dim(im))==4 && dim(im)[4]>1) {
    
    # For each channel
    for(i in seq_len(dim(im)[3])) {
      print(i)
      # z sections as least
      .c <- lapply(seq.int(start, end), function(j) im[, , i, j])
      # Save max in first z section
      im[, , i, 1] <- do.call(FUN, .c)
    }
  }
  
  # Drop 4th dimention
  if(length(dim(im))==4)
    im <- im[, , , 1, drop = TRUE]
  
  # Return first z section
  return(im)
}