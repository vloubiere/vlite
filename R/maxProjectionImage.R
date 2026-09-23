#' Maximum project of an image.
#'
#' @param im Image to process.
#' @param start The first z-section to be included in the max projection.
#' @param end  The last z-section to be included in the max projection.
#'
#' @returns
#' @export
#'
#' @examples
maxProjectionImage <- function(
    im,
    start = 1,
    end = dim(im)[4]
) {
  stopifnot(inherits(im, "AnnotatedImage") | inherits(im, "array"))
  
  # If several z-sections
  if(length(dim(im))==4 && dim(im)[4]>1) {
    
    # For each channel
    for(i in seq_len(dim(im)[3])) {
      print(i)
      # Selected z sections as list
      .c <- lapply(seq.int(start, end), function(j) im[, , i, j])
      # Compute pixel-wise max across selected z and store it in the 1st z slot
      im[, , i, 1] <- do.call(pmax, .c)
    }
  }
  
  # Keep only first 1-section (max) for each color
  if(length(dim(im))==4)
    im <- im[, , , 1, drop = TRUE]
  
  # Return
  return(im)
}