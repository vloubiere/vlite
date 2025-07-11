% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/balloons_plot.R
\name{balloons_plot}
\alias{balloons_plot}
\title{Ballons plot}
\usage{
balloons_plot(
  size.var,
  color.var,
  color.breaks = NULL,
  col = NULL,
  cex = 1,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  gap.axis = 0,
  positive.only = FALSE,
  legend.left.pos = NULL,
  size.legend.title = "Size",
  size.legend.breaks = NULL,
  color.legend.title = "Color"
)
}
\arguments{
\item{size.var}{A color.varrix of numeric values controlling balloons'sizes.}

\item{color.var}{A color.varrix of numeric values controlling balloons'colors.}

\item{color.breaks}{A numeric vector specifying the breakpoints for color mapping.
Defaults to 21 evenly spaced breaks spanning the range of `color.var`.}

\item{col}{A vector of colors corresponding to the values in `breaks`. If `NULL`,
a diverging blue-white-red palette is used for data spanning positive and negative values,
or a sequential viridis-like palette.}

\item{cex}{Expansion factor for balloons. default= 1.}

\item{main}{Title. Default= NA}

\item{gap.axis}{See ?par(). Default= 0}

\item{positive.only}{If set to TRUE, only positive values are shownd. Default= FALSE.}

\item{xaxt}{See ?par(). Default= "s".}

\item{yaxt}{See ?par(). Default= "s".}
}
\value{
Balloon plot
}
\description{
Starting from two color.varrices (dot size and colors, respectively), generates a balloon plot
}
\examples{
color.var <- size.var <- matrix(-3:5, ncol= 3)
rownames(size.var) <- c("r1", "r2", "r3")
colnames(size.var) <- c("A", "B", "C")
balloons_plot(size.var= size.var, color.var= color.var)

}
