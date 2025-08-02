#' Barplot
#'
#' A wrapper around barplot with nicer default settings, allowing input as a list to show both SD and/or individual points.
#'
#' @param x A vector, matrix, data.frame, or a list of numeric variables to plot.
#' @param names.arg Names to plot on the x-axis.
#' @param show.sd If x is a list, should the SD be shown? Default = TRUE.
#' @param show.points If x is a list, should the points be shown? Default = TRUE.
#' @param beside If TRUE and x is a matrix or data.frame, bars are juxtaposed. Default = FALSE (stacked bars).
#' @param xlim X-axis limits. Default = NULL.
#' @param ylim Y-axis limits. Default = NULL.
#' @param col Colors to be used. By default, ?gray.colors will be used.
#' @param sd.arrow.lwd Line width of SD arrows. Default = 0.5.
#' @param sd.arrow.length Length of SD arrows. Defaults to 1/8 of the distance between bar centers.
#' @param pch Plotting character for points. Default = 16.
#' @param pch.col Color for points.
#' @param pch.cex Expansion factor for points.
#' @param pch.jitter Jitter amount for points. Default = 0.2.
#' @param xaxt If tilt.names is set to FALSE, style to be used for the x-axis. Default = "s".
#' @param horiz Not supported at the moment.
#'
#' @return A barplot, and invisibly returns the x and y coordinates of the top of the bars.
#' @export
#'
#' @examples
#' vl_barplot(list(1:3, 5:8), show.sd = FALSE)
#' vl_barplot(matrix(1:6, ncol = 2), beside = TRUE)
vl_barplot <- function(x,
                       names.arg= NULL,
                       show.sd= TRUE,
                       show.points= TRUE,
                       beside= FALSE,
                       width= .6,
                       xlim= NULL,
                       ylim= NULL,
                       col= NULL,
                       sd.arrow.lwd= .5,
                       sd.arrow.length= width/5,
                       pch= 16,
                       pch.col= adjustcolor("grey20", .7),
                       pch.cex= .5,
                       pch.jitter= .2,
                       xaxt= "s",
                       xaxs= "i",
                       frame= F,
                       horiz= FALSE,
                       xlab= NA,
                       ylab= NA,
                       ...)
{
  # Checks ----
  if(horiz)
    stop("horiz not supported atm.")
  if(!is.numeric(c(unlist(x))))
    stop("x should contain numeric values")

  # Simplify list to vector when possible ----
  if(is.list(x) && all(lengths(x)==1))
    x <- unlist(x)

  # Table to matrix ----
  if(is.table(x)) {
    x <- if(is.na(ncol(x)))
      matrix(x, nrow= 1, dimnames = list(names(x), NULL)) else
        matrix(x, ncol= ncol(x), nrow= nrow(x), dimnames = dimnames(x))
  }

  # Data.frame to matrix ----
  if(is.data.frame(x))
    x <- as.matrix(x, drop= F)

  # Compute statistics for lists ----
  if(is.list(x)) {
    # Copy
    x.pts <- x
    # Compute sd
    if(show.sd)
      sd <- sapply(x.pts, sd)
    # Compute mean
    x <- sapply(setNames(x.pts, names(x)), mean)
    # Compute ylim and max var
    if(show.points) {
      yrange <- range(c(0, unlist(x.pts)))
      ymax <- sapply(x.pts, function(x) x[which.max(abs(x))])
    } else if(show.sd) {
      yrange <- range(c(0, x+sd, x-sd))
      ymax <- mapply(function(x, y) ifelse(x>=0, x+y, x-y), x= x, y= sd)
    }
  }

  # Juxtapose columns ----
  if(is.matrix(x) && beside && nrow(x)>1) {
    mat.ncol <- ncol(x)
    x <- c(x)
  } else
    beside <- F

  # Default names.arg ----
  if(is.null(names.arg))
    names.arg <- if(is.matrix(x))
      colnames(x) else
        names(x)

  # Compute at ----
  at <- if(is.matrix(x)) {
    seq(ncol(x))
  } else if(beside) {
    mat.nrow <- length(x)/mat.ncol
    norm.width <- width/mat.nrow
    # Rep per row
    at <- rep(seq(mat.ncol), each= mat.nrow)
    # Space regularly
    at <- at+cumsum(c(0, rep(norm.width, mat.nrow-1)))
    # Shift
    at-width/2+norm.width/2
  } else {
    seq_along(x)
  }

  # Default colors ----
  if(is.null(col)) {
    col <- if(is.matrix(x))
      rev(gray.colors(nrow(x))) else if(beside)
        rev(gray.colors(mat.nrow)) else
          "lightgrey"
  }

  # Compute xlim ----
  if(is.null(xlim)) {
    xlim <- if(is.matrix(x))
      c(0.5, ncol(x)+0.5) else if(beside)
        c(0.5, mat.ncol+0.5) else
          c(0.5, length(x)+0.5)
  }

  # Compute ylim ----
  if(is.null(ylim)) {
    ylim <- if(exists("yrange"))
      yrange else if(is.matrix(x))
        range(c(0, colSums(x))) else
          range(c(0, x))
  }

  # Initiate plot ----
  plot(NA,
       xlim= xlim,
       ylim= ylim,
       xaxt= "n",
       xaxs= xaxs,
       frame= frame,
       xlab= xlab,
       ylab= ylab,
       ...)

  # Plot bars ----
  if(is.matrix(x)) {
    # Stacked barplot
    sapply(seq(ncol(x)), function(i) {
      .c <- cumsum(x[,i])
      rect(i-width/2, c(0, .c[-length(.c)]), i+width/2, .c, col= col)
    })
  } else if(beside) {
    # Juxtaposed barplot
    rect(at-norm.width/2, 0, at+norm.width/2, x, col= col)
    # Only report the position of central bar
    at.center <- seq(mat.ncol)
  } else {
    # Regular barplot
    rect(at-width/2, 0, at+width/2, x, col= col)
  }

  # Add sd ----
  if(show.sd && is.numeric(sd)) {
    segments(at,
             x-sd,
             at,
             x+sd,
             lwd= sd.arrow.lwd,
             xpd= NA)
    segments(at-sd.arrow.length,
             x+sd,
             at+sd.arrow.length,
             x+sd,
             lwd= sd.arrow.lwd,
             xpd= NA)
    segments(at-sd.arrow.length,
             x-sd,
             at+sd.arrow.length,
             x-sd,
             lwd= sd.arrow.lwd,
             xpd= NA)
  }

  # Add points ----
  if(show.points & exists("x.pts")) {
    x.pos <- rep(at, lengths(x.pts))
    y.pos <- unlist(x.pts)
    points(jitter(x.pos, amount = pch.jitter),
           y.pos,
           pch= pch,
           col= unlist(pch.col),
           xpd= NA,
           cex= pch.cex)
  }

  # Add x labels ----
  if(xaxt!="n") {
    check.width <- strwidth(names.arg, units = "user")
    if(any(check.width>width)) {
      tiltAxis(if(beside) at.center else at,
               labels = names.arg)
    } else if(!is.null(names.arg)) {
      axis(1,
           at= if(beside) at.center else at,
           names.arg,
           lwd= 0)
    }
  }

  # Compute max y value if missing ----
  if(!exists("ymax")) {
    ymax <- if(is.matrix(x))
      sapply(colSums(x), function(x) max(c(0, x))) else
        sapply(x, function(x) max(c(0, x)))
  }

  # Return bar positions and top var ----
  invisible(data.table(x= at, y= ymax))
}
