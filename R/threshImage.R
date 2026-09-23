#' Compute automatic intensity thresholds for an image
#'
#' Computes lower and upper intensity cutoffs for image clipping. The lower
#' cutoff is estimated using Otsu's method by default, or optionally from an
#' intensity quantile or the first background density peak. The upper cutoff
#' is quantile-based.
#'
#' @param im Input image or numeric array.
#' @param auto.bg.scaling Scaling factor applied to the upper boundary of the
#' first background intensity peak to define the lower cutoff. A value of
#' 1 uses the boundary directly, while values greater than 1 shift the
#' cutoff toward higher intensities. Default = NULL, corresponding to Otsu thresholding.
#' @param adjust.bandwidth Adjust the bandwidth used to compute density. Higher= smoother. Default= 1
#' @param min.quantile Quantile used as the lower cutoff. Mutually exclusive
#'   with auto.bg.scaling. Default= NULL, corresponding to Otsu thresholding.
#' @param max.quantile Quantile used as the upper cutoff. Default= 0.999.
#'
#' @return A numeric vector containing the lower and upper intensity cutoffs.
#'
#' @export
threshImage <- function(
    im,
    auto.bg.scaling= NULL,
    adjust.bandwidth= 1,
    min.quantile= NULL,
    max.quantile= 0.999
)
{
  # Checks
  if (!is.null(auto.bg.scaling) && !is.null(min.quantile))
    stop("`auto.bg.scaling` and `min.quantile` cannot both be supplied.")
  
  # Compute min clipping value based on otsu/density/quantiles
  min.cutoff <- if(!is.null(min.quantile)) {
    
    # Quantiles method
    quantile(as.numeric(im), min.quantile)
    
  } else if(!is.null(auto.bg.scaling)) {
    
    # Estimate the upper limit of the first (background) density peak
    .c <- as.numeric(c(im))
    d <- stats::density(.c, adjust= adjust.bandwidth, na.rm = TRUE)
    deriv <- diff(d$y) > 0
    r <- data.table::rleidv(deriv)
    d$x[max(which(r == 2))]*auto.bg.scaling
    
    # Method based on background SD
    # bg.limit <- d$x[max(which(r == 2))]
    # # Define the lower cutoff from inferred background pixels
    # bg.pixels <- .c[.c <= bg.limit]
    # mean(bg.pixels, na.rm = TRUE) + auto.bg.sd * sd(bg.pixels, na.rm = TRUE)
    
  } else {
    
    # Defaults to otsu method
    EBImage::otsu(im)
    
  }
  
  # Compute max clipping value based on quantiles
  max.cutoff <- quantile(as.numeric(im), max.quantile)
  
  # Return
  return(c(min.cutoff, max.cutoff))
}
