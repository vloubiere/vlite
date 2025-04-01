#' Plot Matrix as Enhanced Heatmap Image
#'
#' @description
#' A wrapper around base R's image() function that provides enhanced visualization
#' of matrices as heatmaps. The function offers improved default settings and
#' displays the image in a more intuitive orientation for data analysis, with
#' automatic handling of color schemes, NA values, and axis labels.
#'
#' @param mat A numeric matrix to be displayed as a heatmap.
#' @param zlim Numeric vector of length 2 specifying the minimum and maximum values
#'   for color mapping (e.g., c(-1, 1)). If NULL (default), limits are automatically
#'   calculated based on data range. Note: this argument is ignored if breaks are
#'   specified directly.
#' @param breaks Numeric vector of break points for color mapping. Takes precedence
#'   over `zlim` if both are specified. If NULL (default), breaks are automatically
#'   calculated either from `zlim` if provided, or from the data range. Break points
#'   define the boundaries between colors, so if using n colors, breaks should be a
#'   vector of length n+1.
#' @param col Color palette for the heatmap. If NULL (default), automatically
#'   selects between a diverging blue-white-red palette for data (or zlim) centered around
#'   zero, or a sequential blue-yellow palette for non-negative data.
#' @param xlim Numeric vector of length 2 giving x-axis limits. Default is
#'   c(0.5, ncol(mat)+0.5).
#' @param ylim Numeric vector of length 2 giving y-axis limits. Default is
#'   c(nrow(mat)+0.5, 0.5).
#' @param show.rownames Logical; whether to show row names (TRUE by default).
#' @param show.colnames Logical; whether to show column names (TRUE by default).
#' @param tilt.colnames Logical; whether to tilt column names for better
#'   readability (TRUE by default).
#' @param xlab Character string for x-axis label (NA by default).
#' @param ylab Character string for y-axis label (NA by default).
#' @param na.col Color to use for NA values ("darkgrey" by default).
#' @param show.numbers Logical or matrix; if TRUE or a matrix is provided,
#'   displays values in cells. Default is FALSE.
#' @param numbers.cex Numeric scaling factor for number size when show.numbers
#'   is TRUE (0.8 by default).
#'
#' @return A plot of the matrix as a heatmap is produced on the current
#'   graphics device. Invisibly returns the breaks and the colors used.
#'
#' @details
#' The function implements several improvements over base R's image():
#' * Automatic color palette selection based on data characteristics
#' * Intelligent handling of NA values
#' * Better default orientation (rows top-to-bottom)
#' * Optional tilted column labels for better readability
#' * Ability to display numerical values in cells
#' * Automatic handling of row and column names
#'
#' The color scheme is automatically selected based on the data range:
#' * For data spanning positive and negative values: blue-white-red
#' * For non-negative data: blue-yellow
#'
#' @examples
#' # Basic usage with random data
#' set.seed(123)
#' mat <- matrix(rnorm(50), 10, 5)
#' colnames(mat) <- paste0("Sample", 1:5)
#' rownames(mat) <- paste0("Gene", 1:10)
#' vl_image(mat)
#'
#' # Show cell values
#' vl_image(mat, show.numbers = TRUE)
#'
#' # Custom color palette
#' vl_image(mat, col = heat.colors(100))
#'
#' # Non-tilted column names
#' vl_image(mat, tilt.colnames = FALSE)
#'
#' @seealso
#' \code{\link{image}} for the base R image function
#'
#' @export
vl_image <- function(mat,
                   zlim= NULL,
                   breaks= NULL,
                   col= NULL,
                   xlim= c(0.5, ncol(mat)+0.5),
                   ylim= c(nrow(mat)+0.5, 0.5),
                   show.rownames= TRUE,
                   show.colnames= TRUE,
                   tilt.colnames= TRUE,
                   xlab= NA,
                   ylab= NA,
                   legend.cex= 1,
                   legend.title= "Value",
                   na.col= "darkgrey",
                   show.numbers= FALSE,
                   numbers.cex= .8)
{
  # Checks ----
  if(is.null(colnames(mat)))
    colnames(mat) <- seq(ncol(mat))
  if(is.null(rownames(mat)))
    rownames(mat) <- seq(nrow(mat))
  if(isTRUE(show.numbers))
    show.numbers <- mat

  # Default colors ----
  if(is.null(col)) {
    lims <- if(is.null(zlim)) {
      range(mat, na.rm= TRUE)
    } else
      zlim
    col <- if(any(lims>0) && any(lims<0)) {
      # Matrix contains both positive and negative values
      c("royalblue1", "white", "red")
    } else {
      # Or not
      c("blue", "yellow")
    }
    col <- colorRampPalette(col)(ifelse(is.null(breaks), 100, length(breaks)-1))
  }

  # Compute default breaks ----
  if(is.null(breaks)) {
    # If zlim not provided, compute it
    if(is.null(zlim)) {
      zlim <- if(any(mat>0) && any(mat<0)) {
        # Matrix contains both positive and negative values
        lim <- max(abs(mat), na.rm= TRUE)
        c(-lim, lim)
      } else {
        # Or not
        range(mat, na.rm = TRUE)
      }
    }
    # Compute equidistant breaks using zlim
    breaks <- seq(zlim[1], zlim[2], length.out = length(col)+1)
  } else if(!is.null(zlim))
    warning("'zlim' argument is ignored when 'breaks' are specified directly")


  # Transpose matrix ----
  im <- t(mat)

  # Replace outliers and NAs ----
  im[im<min(breaks, na.rm= TRUE)] <- min(breaks, na.rm= TRUE)
  im[im>max(breaks, na.rm= TRUE)] <- max(breaks, na.rm= TRUE)
  im[is.na(im)] <- min(im, na.rm= TRUE)-1

  # Adjust breaks and colors to affor NAs ----
  NAbreaks <- c(breaks[1]-c(2,1), breaks[-1])
  NAcols <- c(na.col, col)

  # Plot image ----
  image(x= seq(nrow(im)),
        y= seq(ncol(im)),
        z= im,
        breaks= NAbreaks,
        col= NAcols,
        xlim= xlim,
        ylim= ylim,
        xlab= xlab,
        ylab= ylab,
        axes= FALSE)

  # Plot axes ----
  if(show.rownames) {
    axis(2,
         seq(nrow(mat)),
         rownames(mat),
         lwd = 0,
         lwd.ticks = 0)
  }
  if(show.colnames) {
    if(tilt.colnames) {
      tiltAxis(x= seq(ncol(mat)),
               labels= colnames(mat))
    } else {
      axis(1,
           seq(ncol(mat)),
           colnames(mat),
           lwd = 0,
           lwd.ticks = 0)
    }
  }

  # Plot numbers ----
  if(is.matrix(show.numbers)) {
    text(x= c(col(show.numbers)),
         y= c(row(show.numbers)),
         labels = c(show.numbers),
         cex= numbers.cex)
  }

  # Return breaks and colors ----
  obj <- list(breaks= breaks, col= col)
  invisible(obj)
}
