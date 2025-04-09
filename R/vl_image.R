#' Plot Matrix as Heatmap
#'
#' A wrapper around base R's `image()` function for enhanced visualization of matrices as heatmaps.
#' The function provides intuitive orientation, automatic handling of color schemes, NA values,
#' and axis labels, along with options for displaying numerical values in cells.
#'
#' @param mat A numeric matrix to be displayed as a heatmap.
#' @param breaks A numeric vector specifying the breakpoints for color mapping.
#'   Defaults to 21 evenly spaced breaks spanning the range of `mat`.
#' @param col A vector of colors corresponding to the values in `breaks`. If `NULL`,
#'   a diverging blue-white-red palette is used for data spanning positive and negative values,
#'   or a sequential blue-yellow palette.
#' @param xlim Numeric vector of length 2 giving x-axis limits. Default is `c(0.5, ncol(mat) + 0.5)`.
#' @param ylim Numeric vector of length 2 giving y-axis limits. Default is `c(nrow(mat) + 0.5, 0.5)`.
#' @param show.rownames Logical; whether to display row names. Default is `TRUE`.
#' @param show.colnames Logical; whether to display column names. Default is `TRUE`.
#' @param tilt.colnames Logical; whether to tilt column names for better readability. Default is `TRUE`.
#' @param xlab Character string for the x-axis label. Default is `NA`.
#' @param ylab Character string for the y-axis label. Default is `NA`.
#' @param na.col Color to use for NA values. Default is `"lightgrey"`.
#' @param show.numbers Logical or matrix; if `TRUE` or a matrix is provided, displays values in cells. Default is `FALSE`.
#' @param numbers.cex Numeric scaling factor for the size of numbers when `show.numbers` is `TRUE`. Default is `0.8`.
#'
#' @return Invisibly returns a list containing the `breaks` and `col` used for the heatmap.
#'
#' @details
#' This function improves upon base R's `image()` by:
#' - Automatically selecting color palettes based on data characteristics.
#' - Handling NA values intelligently.
#' - Displaying row and column names by default.
#' - Allowing optional display of numerical values in cells.
#' - Providing better default orientation (rows top-to-bottom).
#'
#' The default behavior for `breaks` is to create 21 evenly spaced intervals spanning the range of `mat`.
#' If `mat` contains both positive and negative values, the range is symmetrically extended around zero.
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
#' vl_image(mat, col = heat.colors(101))
#'
#' # Non-tilted column names
#' vl_image(mat, tilt.colnames = FALSE)
#'
#' @seealso
#' \code{\link{image}} for the base R `image()` function.
#'
#' @export
vl_image <- function(mat,
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
                     na.col= "lightgrey",
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

  # Default breaks ----
  if(is.null(breaks)) {
    breaks <- if(min(mat, na.rm= TRUE)<0 & max(mat, na.rm= TRUE)>0) {
      lims <- max(abs(mat), na.rm= TRUE)
      # Centered on 0
      seq(-lims,
          lims,
          length.out= 21)
    } else {
      # Otherwise
      seq(min(mat, na.rm= TRUE),
          max(mat, na.rm= TRUE),
          length.out= 21)
    }
  }

  # First interval will contain values <= min break ----
  breaks <- c(breaks[1], breaks)

  # Default colors ----
  if(is.null(col)) {
    col <- if(breaks[1] < 0 & breaks[length(breaks)] > 0){
      # Positive and negative values
      c("royalblue1", "white", "red")
    } else {
      # Otherwise
      c("blue", "yellow")
    }
  }
  col <- colorRampPalette(col)(length(breaks)-1)

  # Transpose matrix ----
  im <- t(mat)

  # Clip outlier values to min/max color breaks ----
  im[im<min(breaks, na.rm= TRUE)] <- min(breaks, na.rm= TRUE)
  im[im>max(breaks, na.rm= TRUE)] <- max(breaks, na.rm= TRUE)

  # Set NA values to min(breaks) -1 ----
  im[is.na(im)] <- min(im, na.rm= TRUE)-1

  # Adjust breaks to afford NAs ----
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
  obj <- list(breaks= breaks[-1],
              col= col)
  invisible(obj)
}
