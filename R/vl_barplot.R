#' Barplot
#'
#' A wrapper around barplot with nicer default settings, allowing input as a list to show both SD and/or individual points.
#'
#' @param x A vector, matrix, data.frame, or a list of numeric variables to plot.
#' @param names.arg Names to plot on the x-axis.
#' @param top.labels Labels to plot on the top of each bar.
#' @param show.sd If x is a list, should the SD be shown? Default = TRUE.
#' @param show.points If x is a list, should the points be shown? Default = TRUE.
#' @param beside If TRUE and x is a matrix or data.frame, bars are juxtaposed. Default = FALSE (stacked bars).
#' @param xlim X-axis limits. Default = NULL.
#' @param ylim Y-axis limits. Default = NULL.
#' @param tilt.names Should bar names be tilted? Default = TRUE.
#' @param top.labels.cex Expansion factor for top.labels.
#' @param sd.arrow.lwd Line width of SD arrows. Default = 0.5.
#' @param sd.arrow.length Length of SD arrows. Defaults to 1/8 of the distance between bar centers.
#' @param pch Plotting character for points. Default = 16.
#' @param pch.col Color for points.
#' @param pch.cex Expansion factor for points.
#' @param pch.jitter Jitter amount for points. Default = 0.2.
#' @param xaxt If tilt.names is set to FALSE, style to be used for the x-axis. Default = "s".
#' @param horiz Not supported at the moment.
#' @param ... Additional arguments passed to \code{barplot()}.
#'
#' @return A barplot.
#' @export
#'
#' @examples
#' vl_barplot(list(1:3, 5:8), top.labels = c("A", "B"), show.sd = FALSE)
#' vl_barplot(matrix(1:6, ncol = 2), top.labels = c("A", "B"), beside = TRUE)
vl_barplot <- function(x,
                       names.arg= NULL,
                       top.labels= NULL,
                       show.sd= TRUE,
                       show.points= TRUE,
                       beside= FALSE,
                       xlim= NULL,
                       ylim= NULL,
                       tilt.names= TRUE,
                       top.labels.cex= .7,
                       sd.arrow.lwd= .5,
                       sd.arrow.length= ifelse(length(x)>1, diff(bar[c(1,2)])/8, .2),
                       pch= 16,
                       pch.col= adjustcolor("grey20", .7),
                       pch.cex= .5,
                       pch.jitter= .2,
                       xaxt= "s",
                       width= .6,
                       space= NULL,
                       xaxs= "i",
                       horiz= FALSE,
                       ...)
{
  # Checks ----
  if(horiz)
    stop("horiz not supported atm.")
  # Default names ----
  if(is.null(names.arg))
    names.arg <- names(x)
  # Table to matrix ----
  if(is.table(x))
    x <- matrix(x, ncol= ncol(x), nrow= nrow(x), dimnames = dimnames(x))

  # If x is a list, compute stats ----
  if(is.list(x)) {
    # Copy
    x.pts <- x
    # Compute sd
    if(show.sd)
      sd <- sapply(x.pts, sd)
    # Compute mean
    x <- sapply(x.pts, mean)
    # Compute ylim and max var
    if(show.points) {
      ylim <- range(c(0, unlist(x.pts)))
      top.var <- sapply(x.pts, function(x) x[which.max(abs(x))])
    } else if(show.sd) {
      ylim <- range(c(0, x+sd, x-sd))
      top.var <- mapply(function(x, y) ifelse(x>=0, x+y, x-y), x= x, y= sd)
    }
  } else {
    # compute top var
    top.var <- if(is.vector(x))
      x else if(beside)
        unlist(c(x)) else
          apply(x, 2, sum, na.rm= T)
  }

  # Default xlim ----
  if(is.null(xlim)) {
    xlim <- if(is.matrix(x)) {
      if(beside & nrow(x)>1) {
        c(0.5, ncol(x)+.5)
      } else
        c(0.5, ncol(x)+.5)
    }else if(is.vector(x))
      c(0.5, length(x)+0.5)
  }

  # Compute ideal space ----
  if(is.null(space)) {
    space <- if(is.matrix(x)) {
      if(beside & nrow(x)>0) {
        width <- width/nrow(x)
        message("beside= TRUE -> bar width adjusted")
        .c <- c((1-width/2)/width, rep((1-(width*nrow(x)))/width, ncol(x)-1))
        c(sapply(.c, function(y) c(y, rep(0, nrow(x)-1))))
      } else
        c((1-width/2)/width, rep((1-width)/width, ncol(x)-1))
    }else if(is.vector(x))
      c((1-width/2)/width, rep((1-width)/width, length(x)-1))
  }

  # Plot ----
  bar <- barplot(x,
                 xlim= xlim,
                 ylim= ylim,
                 beside= beside,
                 xaxs= xaxs,
                 xaxt= ifelse(tilt.names, "n", xaxt),
                 names.arg= names.arg,
                 width= width,
                 space= space,
                 ...)

  # Add x labels ----
  if(tilt.names & xaxt!="n")
    tiltAxis(bar,
             labels = names.arg)

  # Add sd ----
  if(show.sd && is.numeric(sd)) {
    segments(bar,
             x-sd,
             bar,
             x+sd,
             lwd= sd.arrow.lwd,
             xpd= NA)
    segments(bar-sd.arrow.length,
             x+sd,
             bar+sd.arrow.length,
             x+sd,
             lwd= sd.arrow.lwd,
             xpd= NA)
    segments(bar-sd.arrow.length,
             x-sd,
             bar+sd.arrow.length,
             x-sd,
             lwd= sd.arrow.lwd,
             xpd= NA)
  }

  # Add points ----
  if(show.points & exists("x.pts")) {
    x <- rep(bar, lengths(x.pts))
    y <- unlist(x.pts)
    points(jitter(x, amount = pch.jitter),
           y,
           pch= pch,
           col= unlist(pch.col),
           xpd= NA,
           cex= pch.cex)
  }

  # Add bar labels ----
  if(!is.null(top.labels))
  {
    text(bar,
         top.var,
         top.labels,
         pos= ifelse(top.var>=0, 3, 1),
         cex= top.labels.cex,
         xpd= NA)
  }
}
