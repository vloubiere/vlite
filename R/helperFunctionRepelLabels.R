
# Function that takes as input rectangle coordinates and segment coordinates
# and returns the coordinates of a clipped, non-overlapping segment
find_closest_rect_segment <- function(xleft,
                                      ybottom,
                                      xright,
                                      ytop,
                                      seg.x0,
                                      seg.y0,
                                      seg.x1,
                                      seg.y1,
                                      n_points = 10) # Number of points to be used. The higher the more precise the clipping will be
{
  # Helper for between
  between <- function(x, left, right) x >= left & x <= right

  # Compute points along segment
  seg.x <- seq(seg.x0, seg.x1, length.out = n_points)
  seg.y <- seq(seg.y0, seg.y1, length.out = n_points)

  # Find which endpoint is outside the rectangle
  if (between(seg.x0, xleft, xright) & between(seg.y0, ybottom, ytop)) {
    outside.x <- seg.x1
    outside.y <- seg.y1
  } else {
    outside.x <- seg.x0
    outside.y <- seg.y0
  }

  # Compute candidate points on rectangle border
  x.rect <- numeric()
  y.rect <- numeric()
  if (!between(outside.y, ybottom, ytop)) {
    x.rect <- seq(xleft, xright, length.out = n_points)
    y.rect <- if (outside.y <= ybottom) rep(ybottom, n_points) else rep(ytop, n_points)
  }
  if (!between(outside.x, xleft, xright)) {
    y.rect <- c(y.rect, seq(ybottom, ytop, length.out = n_points))
    x.rect <- c(x.rect, if (outside.x <= xleft) rep(xleft, n_points) else rep(xright, n_points))
  }

  # Compute pairwise distances
  min_dist_per_rect <- sapply(seq_along(x.rect), function(i) {
    dists <- sqrt((x.rect[i] - seg.x)^2 + (y.rect[i] - seg.y)^2)
    min(dists)
  })

  # Find the index of the closest rect point
  closest <- which.min(min_dist_per_rect)

  # Return the coordinates of the segment of interest
  list(
    x0 = outside.x,
    y0 = outside.y,
    x1 = x.rect[closest],
    y1 = y.rect[closest]
  )
}
