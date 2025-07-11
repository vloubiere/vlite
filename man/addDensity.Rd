% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addDensity.R
\name{addDensity}
\alias{addDensity}
\title{Add Marginal Density Plots to a Scatterplot}
\usage{
addDensity(
  x,
  y = NULL,
  col = adjustcolor("grey", 0.6),
  group = col,
  border = NA,
  height = 1,
  width = 1
)
}
\arguments{
\item{x}{Numeric vector for the x-axis.}

\item{y}{Numeric vector for the y-axis. If NULL, y is set to x and x becomes the index.}

\item{col}{Fill color(s) for the density polygons. Must be length 1 or match the length of x. Default= adjustcolor("grey", 0.6).}

\item{group}{Optional grouping factor. Used to sort data before plotting. Must be length 1 or match x. Default= col argument.}

\item{border}{Border color(s) for the polygons. Must be length 1 or match x. Default is NA.}

\item{height}{Height of the x-axis density (in lines).}

\item{width}{Width of the y-axis density (in lines).}
}
\value{
Adds density polygons to the current plot. No return value.
}
\description{
Adds density polygons to the top and right sides of an existing scatterplot, optionally grouped and colored.
}
\examples{
set.seed(1)
dat <- data.table(x= rnorm(1000, 1, 0.5),
                  y= rnorm(1000, 1, 0.5))
setorder(dat, "x")
dat[, col:= rep(c("blue", "red"), each = 500)]
vl_par()
plot(dat$x, dat$y, col = dat$col)
addDensity(x = dat$x, y = dat$y, col = adjustcolor(dat$col, 0.3))
}
