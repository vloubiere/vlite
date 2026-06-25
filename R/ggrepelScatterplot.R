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
    label.col = "black",
    xlim= NULL,
    ylim= NULL,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    cex = 1,
    cex.pts= 1,
    cex.lab = 1,
    plot= T
) {
  # Scaling factor
  cex.pts <- cex.pts*1.5
  cex.lab <- cex.lab*3
  
  # Checks
  stopifnot(length(x) == length(y),
            length(x) == length(col),
            length(x) == length(label))
  
  # Default labels and colors
  if (is.null(label.col)) label.col <- col
  if (length(label.col) == 1L) label.col <- rep(label.col, length(x))
  stopifnot(length(label.col) == length(x))
  
  # Default limits
  if(is.null(xlim)) {
    xlim <- range(x, na.rm = T)
    exp <- diff(xlim)*0.04
    xlim[1] <- xlim[1]-exp
    xlim[2] <- xlim[2]+exp
  }
  if(is.null(ylim)) {
    ylim <- range(y, na.rm = T)
    exp <- diff(ylim)*0.04
    ylim[1] <- ylim[1]-exp
    ylim[2] <- ylim[2]+exp
  }
  
  # Make DF
  df <- data.frame(
    x = x,
    y = y,
    col = col,
    label = label,
    label.col = label.col,
    stringsAsFactors = FALSE
  )
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
  })
  
  # Compute plot
  pl <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = col), size = cex.pts * cex) +
    scale_color_identity() +
    geom_text_repel(aes(label = label, color = label.col),
                    size = cex.lab * cex,
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    min.segment.length = 0.1) +
    labs(x = xlab, y = ylab) +
    coord_cartesian(
      xlim = xlim,
      ylim = ylim,
      clip = "off"
    ) +
    ggtitle(main) +
    theme_classic()
  
  # Plot
  if(plot)
    plot(pl) else
      return(pl)
}