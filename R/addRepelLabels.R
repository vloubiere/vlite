#' addRepelLabels
#'
#' Add repelled labels to an existing scatterplot.
#'
#' @param x x coordinates.
#' @param y y coordinates.
#' @param labels labels to plot.
#' @param label.col A vector of colors for labels. Default= "black".
#' @param label.cex Expansion factor for labels.
#' @param seg.col A vector of colors for segments. Default= "black".
#' @param rect.draw Draw rectangle boxes around labels? Default= TRUE.
#' @param rect.col A vector of colors for boxes. Default= adjustcolor("white", .3).
#' @param no.overlap If set to TRUE, overlapping names will be removed. Default= FALSE.
#' @param point_size Approximate size of points, defined as a fraction of char size.
#' Used in conjunction with 'point_padding' x and y (see below). Default= 0.5
#' @param point_padding_x Extra margin around data points, defined as a fraction of char size. Default= 0.
#' @param point_padding_y Extra margin around data points, defined as a fraction of char size. Default= 0.
#' @param box_padding_x Extra margin around text bounding boxes, defined as a fraction of char size. Default= 0.5.
#' @param box_padding_y Extra margin around text bounding boxes, defined as a fraction of char size. Default= 0.5.
#' @param force_push Magnitude of the push force. Default= 1e-5.
#' @param force_pull Magnitude of the pull force. Default= 1e-5.
#'
#' @return
#' @export
#'
#' @examples
#' data("mtcars")
#' mtcars$car <- rownames(mtcars)
#'
#' plot(x = mtcars$wt, y = mtcars$mpg)
#'
#' addRepelLabels(x = mtcars$wt,
#' y = mtcars$mpg,
#' labels = mtcars$car,
#' cex = 0.7)
#'
addRepelLabels <- function(x,
                           y= NULL,
                           labels,
                           cex= 0.7,
                           label.col= "black",
                           seg.col= "black",
                           rect.draw= TRUE,
                           rect.col= adjustcolor("white", .3),
                           no.overlap= FALSE,
                           point_size= 0.5,
                           points_padding_x = 0,
                           points_padding_y = 0,
                           box_padding_x = .5,
                           box_padding_y = .5,
                           force_push = 1e-05,
                           force_pull = 1e-05)
{
  # Checks ----
  if(is.null(y))
  {
    y <- x
    x <- seq(y)
  }

  # Initiate table ----
  dat <- data.table(x,
                    y,
                    labels,
                    col= label.col,
                    cex= cex,
                    rect.col= rect.col,
                    seg.col= seg.col)
  dat <- dat[!is.na(x) & !is.na(y)]
  dat[, width:= strwidth(labels, cex= cex)]
  dat[, height:= strheight(labels, cex= cex)]
  dat[, x1:= x-width/2-strwidth("M")*box_padding_x*cex]
  dat[, y1:= y-height/2-strheight("M")*box_padding_y*cex]
  dat[, x2:= x+width/2+strwidth("M")*box_padding_x*cex]
  dat[, y2:= y+height/2+strheight("M")*box_padding_y*cex]

  # Points and boxes DF ----
  posdf <- as.data.frame(dat[, .(x, y)])
  boxdf <- as.data.frame(dat[, x1:y2])

  # Compute positions boxes
  moved <- latticetools::repel_boxes(data_points = as.matrix(posdf),
                                     point_size = rep(point_size, nrow(posdf)),# npc unit
                                     point_padding_x = points_padding_x,
                                     point_padding_y = points_padding_y,
                                     boxes = as.matrix(boxdf),
                                     xlim = par("usr")[c(1,2)],
                                     ylim = par("usr")[c(3,4)],
                                     hjust = rep(0, nrow(posdf)),
                                     vjust = rep(0, nrow(posdf)),
                                     force_push = force_pull,
                                     force_pull = force_pull,
                                     max_time = 3,
                                     max_overlaps = 10,
                                     max_iter = 10000,
                                     direction = "both",
                                     verbose = FALSE)
  setnames(moved, c("x", "y"), c("x_moved", "y_moved"))
  dat <- cbind(dat, moved)
  dat[, x1:= x1+(x_moved-x)]
  dat[, x2:= x2+(x_moved-x)]
  dat[, y1:= y1+(y_moved-y)]
  dat[, y2:= y2+(y_moved-y)]

  # Remove overlapping boxes
  if(no.overlap && any(dat$too_many_overlaps))
  {
    removed <- sum(dat$too_many_overlaps, na.rm= T)
    warning(paste0(removed, " overlapping labels were removed from the plot!"))
    dat <- dat[!(too_many_overlaps)]
  }

  # Compute clipped segments ----
  dat[, c("seg.x0", "seg.y0", "seg.x1", "seg.y1"):= {
    find_closest_rect_segment(xleft = x1[1],
                              ybottom = y1[1],
                              xright = x2[1],
                              ytop = y2[1],
                              seg.x0 = x[1],
                              seg.y0 = y[1],
                              seg.x1 = x_moved[1],
                              seg.y1 = y_moved[1],
                              n_points = 10)
  }, .(x1, y1, x2, y2, x, y, x_moved, y_moved)]

  # # Check if segments are long enough ----
  # dat[, seg_length.x:= min(abs(c(seg.x0, seg.x1))), .(seg.x0, seg.x1)] # Check if segments are long enough
  # dat[, seg_length.y:= min(abs(c(seg.y0, seg.y1))), .(seg.y0, seg.y1)]
  # dat[, seg_length:= seg_length.x>strwidth("M", cex= cex/4) | seg_length.y>strheight("M", cex= cex/4)]

  # Plot segments ----
  with(dat, segments(seg.x0, seg.y0, seg.x1, seg.y1, col = seg.col))

  # Plot rectangles ----
  if(rect.draw)
    with(dat, rect(x1, y1, x2, y2, col = rect.col, border= NA))

  # Plot labels ----
  with(dat, text(x_moved, y_moved, labels = labels, col = label.col, cex = cex))

  # Return object
  invisible(dat)
}
