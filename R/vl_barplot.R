#' Title
#'
#' @param x
#' @param names.arg
#' @param show.sd
#' @param show.points
#' @param beside
#' @param width
#' @param xlim
#' @param ylim
#' @param col
#' @param sd.arrow.lwd
#' @param sd.arrow.length
#' @param pch
#' @param pch.col
#' @param pch.cex
#' @param pch.jitter
#' @param xaxt
#' @param xaxs
#' @param frame
#' @param horiz
#' @param xlab
#' @param ylab
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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
                       add= F,
                       xlab= "",
                       ylab= "",
                       ...)
{
  # Checks ----
  if(horiz)
    stop("horiz not supported atm.")
  if(!is.numeric(c(unlist(x))))
    stop("x should contain numeric values")

  # Initialize variables to avoid scope issues
  sd <- NULL
  x.pts <- NULL
  mat.nrow <- NULL
  mat.ncol <- NULL
  ymax <- NULL

  # Simplify list to vector when possible ----
  if(is.list(x) && all(lengths(x)==1))
    x <- unlist(x)

  # Table to vector or matrix ----
  if(is.table(x)) {
    x <- if(is.na(ncol(x))) {
      c(x)
    } else
      matrix(x, ncol= ncol(x), nrow= nrow(x), dimnames = dimnames(x))
  }

  # Data.frame to matrix ----
  if(is.data.frame(x))
    x <- as.matrix(x, drop= F)

  # Compute statistics for lists ----
  if(is.list(x)) {
    x.pts <- x
    if(show.sd)
      sd <- sapply(x.pts, sd)
    x <- sapply(setNames(x.pts, names(x)), mean)
    if(show.points) {
      if(is.null(ylim))
        ylim <- range(c(0, unlist(x.pts)))
      ymax <- sapply(x.pts, function(x) x[which.max(abs(x))])
    } else if(!is.null(sd)) {
      if(is.null(ylim))
        ylim <- range(c(0, x+sd, x-sd))
      ymax <- mapply(function(x, y) ifelse(x>=0, x+y, x-y), x= x, y= sd)
    }
  }

  # Juxtapose columns ----
  if(is.matrix(x) && beside && nrow(x)>1) {
    mat.ncol <- ncol(x)
    mat.nrow <- nrow(x)
    x <- c(x)
  } else
    beside <- F

  # At this point, x will remain as a matrix only for stacked barplots
  # In all other cases, x is a vector

  # Default names.arg ----
  if(is.null(names.arg))
    names.arg <- if(is.matrix(x))
      colnames(x) else
        names(x)

  # Compute at ----
  at <- if(is.matrix(x)) {
    seq(ncol(x))
  } else if(beside) {
    norm.width <- width/mat.nrow
    at <- rep(seq(mat.ncol), each= mat.nrow)
    at <- at+cumsum(c(0, rep(norm.width, mat.nrow-1)))
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
    ylim <- if(is.matrix(x))
      range(c(0, colSums(x))) else
        range(c(0, x))
  }

  # Initiate plot ----
  if(isFALSE(add)) {
    plot(NA,
         xlim= xlim,
         ylim= ylim,
         xaxt= "n",
         xaxs= xaxs,
         frame= frame,
         xlab= xlab,
         ylab= ylab,
         ...)
  }

  # Plot bars ----
  if(is.matrix(x)) {
    # Stacked barplot
    cs <- apply(x, 2, function(x) cumsum(c(0, x)))
    x.b <- rep(seq(ncol(x)), each= nrow(x))
    rect(x.b-width/2, unlist(cs[-1,]), x.b+width/2, unlist(cs[-nrow(cs),]), col= col)
  } else if(beside) {
    # Juxtaposed bars
    rect(at-norm.width/2, 0, at+norm.width/2, x, col= col)
    at.center <- seq(mat.ncol)
  } else {
    # Regular barplot
    rect(at-width/2, 0, at+width/2, x, col= col)
  }

  # Add sd ----
  if(!is.null(sd) && is.numeric(sd)) {
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
  if(show.points & !is.null(x.pts)) {
    x.pos <- rep(at, lengths(x.pts))
    y.pos <- unlist(x.pts)
    # Ensure pch.col is a vector of correct length
    if(length(pch.col) == 1) pch.col <- rep(pch.col, length(y.pos))
    points(jitter(x.pos, amount = pch.jitter),
           y.pos,
           pch= pch,
           col= pch.col,
           xpd= NA,
           cex= pch.cex)
  }

  # Add x labels ----
  if(xaxt!="n") {
    check.width <- strwidth(names.arg, units = "user")
    if(any(check.width>width)) {
      if(exists("tiltAxis", mode="function")) {
        tiltAxis(if(beside) at.center else at,
                 labels = names.arg)
      } else {
        warning("tiltAxis function not found; cannot tilt axis labels.")
        axis(1,
             at= if(beside) at.center else at,
             names.arg,
             lwd= 0)
      }
    } else if(!is.null(names.arg)) {
      axis(1,
           at= if(beside) at.center else at,
           names.arg,
           lwd= 0)
    }
  }

  # Compute max y value if missing ----
  if(is.null(ymax)) {
    ymax <- if(is.matrix(x))
      sapply(colSums(x), function(x) max(c(0, x))) else
        sapply(x, function(x) max(c(0, x)))
  }

  # Return bar positions and top var ----
  invisible(data.frame(x= at, y= ymax))
}
