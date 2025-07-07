#' Add Marginal Density Plots to a Scatterplot
#'
#' Adds density polygons to the top and right sides of an existing scatterplot, optionally grouped and colored.
#'
#' @param x Numeric vector for the x-axis.
#' @param y Numeric vector for the y-axis. If NULL, y is set to x and x becomes the index.
#' @param col Fill color(s) for the density polygons. Must be length 1 or match the length of x. Default= adjustcolor("grey", 0.6).
#' @param group Optional grouping factor. Used to sort data before plotting. Must be length 1 or match x. Default= col argument.
#' @param border Border color(s) for the polygons. Must be length 1 or match x. Default is NA.
#' @param height Height of the x-axis density (in lines).
#' @param width Width of the y-axis density (in lines).
#'
#' @return Adds density polygons to the current plot. No return value.
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- data.table(x= rnorm(1000, 1, 0.5),
#'                   y= rnorm(1000, 1, 0.5))
#' setorder(dat, "x")
#' dat[, col:= rep(c("blue", "red"), each = 500)]
#' vl_par()
#' plot(dat$x, dat$y, col = dat$col)
#' addDensity(x = dat$x, y = dat$y, col = adjustcolor(dat$col, 0.3))
addDensity <- function(x,
                       y= NULL,
                       col= adjustcolor("grey", .6),
                       group= col,
                       border= NA,
                       height= 1,
                       width= 1)
{
  # Checks ----
  if(is.null(y)) {
    y <- x
    x <- seq_along(y)
  }

  # Make data table ----
  dat <- data.table(x, y, group, col, border)

  # Compute thickness ----
  height <- diff(grconvertY(c(0,height), "line", "user"))
  width <- diff(grconvertX(c(0,width), "line", "user"))

  # Order ----
  setorderv(dat, "group")

  # Plot ----
  dat[, {
    # X axis
    .d <- density(x = x,
                  from= par("usr")[1],
                  to= par("usr")[2],
                  na.rm= TRUE)
    .d$y <- .d$y-min(.d$y)
    .d$y <- .d$y/max(.d$y)*height+par("usr")[4]
    polygon(x = c(par("usr")[1], .d$x, par("usr")[2]),
            y = c(par("usr")[4], .d$y, par("usr")[4]),
            xpd= NA,
            col= col[1],
            border= border[1])
    # Y axis
    .d <- density(x = y,
                  from= par("usr")[3],
                  to= par("usr")[4],
                  na.rm= TRUE)
    .d$y <- .d$y-min(.d$y)
    .d$y <- .d$y/max(.d$y)*height+par("usr")[2]
    polygon(x = c(par("usr")[2], .d$y, par("usr")[2]),
            y = c(par("usr")[3], .d$x, par("usr")[4]),
            xpd= NA,
            col= col[1],
            border= border[1])
  }, .(group, col)]

}
