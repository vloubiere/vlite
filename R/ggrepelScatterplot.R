#' ggrepelScatterplot
#' 
#' A wrapper around ggrepel to plot a scatterplot with repel labels, but with parameters
#' matching the names of default R graphics.
#'
#' @param x x variable.
#' @param y y variable.
#' @param col Color of the dots.
#' @param label Labels to be added.
#' @param label.col Label color.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param main title.
#' @param cex Expansion factor for point size.
#' @param point_size point_size parameter ggrepel.
#' @param label_size label_size parameter ggrepel.
#'
#' @returns
#' @export
#'
#' @examples
ggrepelScatterplot <- function(
    x,
    y,
    col,
    label = rep("", length(x)),
    label.col = NULL,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    cex = 1,
    point_size = 2.5,
    label_size = 3
) {
  stopifnot(length(x) == length(y),
            length(x) == length(col),
            length(x) == length(label))
  
  if (is.null(label.col)) label.col <- col
  if (length(label.col) == 1L) label.col <- rep(label.col, length(x))
  stopifnot(length(label.col) == length(x))
  
  df <- data.frame(x = x, y = y, col = col, label = label, label.col = label.col,
                   stringsAsFactors = FALSE)
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
  })
  
  ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = col), size = point_size * cex) +
    scale_color_identity() +
    geom_text_repel(aes(label = label, color = label.col),
                    size = label_size, show.legend = FALSE) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) +
    theme_classic()
}