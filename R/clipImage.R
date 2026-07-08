#' Title
#'
#' @param im Image to clip.
#' @param channels Which channels should be clipped (by default, they all will).
#' @param min.bg.sd The number of sd to be added to the mean value of the background peak 
#' to define the minimum threshold. Default= 2.
#' @param min.quantile If specificied, minimum cutoff will be computed based on quantiles. Default= NA.
#' @param max.quantile The maximum quantile to which values will be clipped. Default= 0.99.
#' @param min.value If specified, hard thresholding will be performed (instead of quantile based). Default= NULL.
#' @param max.value If specified, hard thresholding will be performed (instead of quantile based). Default= NULL.
#'
#' @returns
#' @export
#'
#' @examples
clipImage <- function(
    im,
    channels= NULL,
    min.bg.sd= 2,
    min.quantile= NA,
    max.quantile= 0.99,
    min.value= NULL,
    max.value= NULL
) {
  stopifnot(inherits(im, "AnnotatedImage") | inherits(im, "array"))
  
  # Select channels to process
  if(is.null(channels))
    channels <- seq_len(dim(im)[3])
  stopifnot(all(channels <= dim(im)[3]))
  
  # For each channel
  for(i in channels) {
    .c <- if(length(dim(im)) == 4)
      im[, , i, ] else
        im[, , i]
    
    # Compute min clipping value based on density
    min.cutoff <- if(is.null(min.value)) {
      if(is.na(min.quantile)) {
        # Density
        d <- density(c(.c)) 
        # Extract values from the first (bg) peak
        deriv <- sign(diff(d$y)) == 1 # Derivative
        r <- data.table::rleidv(deriv) # rle
        bg <- d$x[1:max(which(r==2))] # bg values (within 1st peak)
        # Compute threshold
        mean(bg)+sd(bg)*min.bg.sd
      } else {
        quantile(as.numeric(.c), min.quantile)
      }
    } else
      min.value

    # Compute max clipping values
    max.cutoff <- if(is.null(max.value))
      quantile(as.numeric(.c), max.quantile) else
        max.value
    
    # Print limits
    print(
      paste0(
        "Clipping value for channel ", i, ": ",
        formatC(min.cutoff, format = "e", digits = 1),
        " / ",
        formatC(max.cutoff, format = "e", digits = 1)
      )
    )
    
    # Clip
    if(max.cutoff>min.cutoff) {
      .c <- (.c - min.cutoff) / diff(c(min.cutoff, max.cutoff))
      .c[.c < 0] <- 0
      .c[.c > 1] <- 1
    }
    
    # Replace
    if(length(dim(im)) == 4)
      im[, , i, ] <- .c else
        im[, , i] <- .c
  }
  
  return(im)
}