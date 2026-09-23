#' Clip and rescale image intensities
#'
#' Clips selected image channels between lower and upper intensity cutoffs,
#' then linearly rescales the retained range to [0, 1]. Cutoffs can be
#' supplied directly or estimated automatically using ?threshImage.
#'
#' @param im Input image or numeric array, with channels stored along the third
#'   dimension and, optionally, z-sections along the fourth dimension.
#' @param channels Integer vector specifying the channels to process. By
#'   default, all channels are processed.
#' @param auto.bg.scaling Scaling factor applied to the upper boundary of the
#' first background intensity peak to define the lower cutoff. A value of
#' 1 uses the boundary directly, while values greater than 1 shift the
#' cutoff toward higher intensities. Default= NULL, corresponding to Otsu thresholding.
#' @param adjust.bandwidth Adjust the bandwidth used to compute density. Higher= smoother. Default= 1
#' @param min.quantile Quantile used as the automatic lower cutoff. Mutually
#'   exclusive with auto.bg.scaling. Default= NULL, corresponding to Otsu
#'   thresholding.
#' @param max.quantile Quantile used as the automatic upper cutoff. Default= 0.999.
#' @param min.value Optional fixed lower cutoff, overriding its automatic
#'   estimate. Default= NULL.
#' @param max.value Optional fixed upper cutoff, overriding its automatic
#'   estimate. Default= NULL.
#'
#' @return The image with selected channels clipped and rescaled to [0, 1].
#'
#' @export
clipImage <- function(
    im,
    channels= NULL,
    auto.bg.scaling= NULL,
    adjust.bandwidth= 1,
    min.quantile= NULL,
    max.quantile= 0.999,
    min.value= NULL,
    max.value= NULL
) {
  stopifnot(inherits(im, "AnnotatedImage") | inherits(im, "array"))
  stopifnot(all(channels >= 1 & channels <= dim(im)[3]))
  
  # Select channels to process
  if(is.null(channels))
    channels <- seq_len(dim(im)[3])
  stopifnot(all(channels <= dim(im)[3]))
  
  # For each channel
  for(i in channels) {
    .c <- if(length(dim(im)) == 4)
      im[, , i, ] else
        im[, , i]
    
    # Compute automatic clipping values
    if(is.null(min.value) | is.null(max.value)) {
      auto.thresh <- threshImage(
        .c,
        auto.bg.scaling = auto.bg.scaling,
        adjust.bandwidth= adjust.bandwidth,
        min.quantile= min.quantile,
        max.quantile= max.quantile
      )
    }
    
    # Retrieve clipping values
    min.cutoff <- if(is.null(min.value))
      auto.thresh[1] else
        min.value
    max.cutoff <- if(is.null(max.value))
      auto.thresh[2] else
        max.value
    
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