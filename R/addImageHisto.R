#' Add intensity-density plots to an image
#'
#' Adds intensity-density plots for selected image channels to the active
#' graphics device. The intensity range retained after thresholding is
#' highlighted for each channel.
#'
#' @param im Input image or numeric array, with channels stored along the third
#'   dimension and, optionally, z-sections along the fourth dimension.
#' @param channels Integer vector specifying the channels to display. By
#'   default, all channels are displayed.
#' @param pos Numeric vector specifying the vertical position of each density
#'   plot. By default, channels are displayed successively from top to bottom.
#' @param col Color used to highlight intensities between the lower and upper
#'   cutoffs. A single color or one color per channel can be supplied.
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
#' @param ymax Optional upper density limit shared across channels. By default,
#'   it is derived from the density above the lower cutoff.
#' @param ... Extra parameters to be passed to lines().
#'
#' @return
#'
#' @export
addImageHisto <- function(
    im,
    channels= NULL,
    pos= NULL,
    col= "red",
    auto.bg.scaling= NULL,
    adjust.bandwidth= 1,
    min.quantile= NULL,
    max.quantile= 0.999,
    min.value= NULL,
    max.value= NULL,
    ymax= NULL,
    ...
)
{
  stopifnot(inherits(im, "AnnotatedImage") | inherits(im, "array"))
  
  # Select channels to process
  if(is.null(channels))
    channels <- seq_len(dim(im)[3])
  stopifnot(all(channels <= dim(im)[3]))
  
  # Plotting positions
  if(is.null(pos))
    pos <- seq_along(channels)
  
  # Colors
  if(length(col)==1 & length(channels)>1)
    col <- rep(col, length(channels))
  
  # Compute fixed param
  width <- diff(par("usr")[c(1,2)])/7
  left <- par("usr")[2]-width*1.1
  height <- diff(par("usr")[c(3,4)])/12
  
  # For each channel
  for(idx in seq_along(channels)) {
    .c <- if(length(dim(im)) == 4)
      im[, , channels[idx], ] else
        im[, , channels[idx]]
    
    # Compute automatic clipping values
    if(is.null(min.value) | is.null(max.value))
    {
      auto.thresh <- threshImage(
        .c,
        auto.bg.scaling= auto.bg.scaling,
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
    
    # Compute density
    d <- density(
      .c, 
      from = 0,
      to = max.cutoff,
      na.rm= T,
      adjust= adjust.bandwidth
    )
    
    # Compute ylim value based on min.cutoff
    clip.y <- if(is.null(ymax)) {
      if(any(d$x>=min.cutoff))
        max(d$y[d$x>=min.cutoff], na.rm= T)*1.5 else
          max(d$y)
    } else
      ymax
    d$y[d$y>clip.y] <- clip.y
    
    # Plotting parameters
    x <- d$x/max(d$x)
    y <- d$y/max(d$y)
    bottom <- par("usr")[4]-height*pos[idx]*1.1
    sel <- (d$x>=min.cutoff & d$x<=max.cutoff) # ploted intensity range
    
    # Plot
    # plot(c(1, 2048), c(2048,1)) # For tests
    if(any(sel))
      segments(
        x0 = left+x[sel]*width,
        x1 = left+x[sel]*width,
        y0 = bottom,
        y1 = bottom+y[sel]*height,
        col= col[idx]
      )
    lines(left+x*width, bottom+y*height, col= "white", ...)
  }
}