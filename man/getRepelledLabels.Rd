% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repelScatterplot.R
\name{getRepelledLabels}
\alias{getRepelledLabels}
\title{Get bounding boxes for repelled text labels}
\usage{
getRepelledLabels(
  x,
  y,
  labels,
  cex,
  xlim,
  ylim,
  point_size,
  box_padding_x,
  box_padding_y,
  points_padding_x,
  points_padding_y,
  force_push,
  force_pull
)
}
\arguments{
\item{x, y}{(numeric) vectors representing x and y coordinates of points}

\item{labels}{(character) vector of labels}

\item{cex}{(numeric) size of text labels}

\item{xlim}{numeric vector representing the limits on the x axis like
\code{c(xmin, xmax)}}

\item{ylim}{numeric vector representing the limits on the y axis like
 \code{c(ymin, ymax)}
‘force_push’ ‘force_pull’ ‘transform_coord’}

\item{point_size}{(numeric) approximate size of points in npc relative cooordinates. Can be small,
is used in conjunction with "point_padding" options}

\item{box_padding_x}{(numeric) extra margin around text bounding boxes}

\item{box_padding_y}{(numeric) extra margin around text bounding boxes}

\item{force_push}{(numeric) magnitude of the push force}

\item{force_pull}{(numeric) magnitude of the pull force}

\item{point_padding_x}{(numeric) extra margin around data points}

\item{point_padding_y}{(numeric) extra margin around data points}

\item{transform_coord}{(character) the coordinate system used for calculating
bounding boxes. The default, "npc", is to transfrom all cordinates to
'Normalised Parent Coordinates'. When calling externally, use "native"
instead, see examples.}

\item{...}{other arguments passed to the function}
}
\description{
Get bounding boxes for repelled text labels
}
\examples{
data("mtcars")
mtcars$car <- rownames(mtcars)

# with base graphics
# ---------------------------
plot(mtcars$mpg ~ mtcars$wt, pch = 19)

coords <- getRepelledLabels(
  x = mtcars$wt, y = mtcars$mpg,
  labels = mtcars$car,
  cex = 0.7, point_size = 0.01,
  xlim = range(mtcars$wt),
  ylim = range(mtcars$mpg),
  box_padding_x = 0.05,
  box_padding_y = 0.05,
  point_padding_x = 0.01,
  point_padding_y = 0.2,
  force_push = 1e-05,
  force_pull = 1e-05,
  transform_coord = "native"
)

with(coords, rect(x1_box, y1_box, x2_box, y2_box, col = "white"))
with(coords, text(x, y, labels = label, cex = 0.7))

# with ggplot2
# ---------------------------

\dontrun{
library(ggplot2)
ggplot(mtcars, aes(y = mpg, x = wt)) +
  geom_point() +
  geom_label(aes(x = coords$x, y = coords$y, label = coords$label))
}

}
