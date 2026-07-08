#' Title
#'
#' @param channels A vector of names for the channels (DAPI, GFP...).
#' @param col Channel colors
#' @param title A title for the image.
#' @param cex Expansion factor applied to test.
#'
#' @returns
#' @export
#'
#' @examples
addImageLegend <- function(channels= "DAPI", col= "blue", title, cex= 1) {
  
  # Adjust cex factor
  cex <- cex*2
  
  # Compute starting x and y
  x <- par("usr")[1]
  y <- par("usr")[4]-strheight("M")*cex
  
  # Plot title first (in white)
  if(!missing(title)) {
    text(x, y, title, col= "white", pos= 4, cex= cex)
    y <- y-strheight("M")*cex*1.5
  }
  
  # Plot channel names
  for(i in seq(length(channels))) {
    text(x, y, channels[i], col= col[i], pos= 4, cex= cex)
    y <- y-strheight("M")*cex*1.5
  }
}