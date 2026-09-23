#' vl_barplot
#' 
#' A wrapper around the barplot function where the bars are regularly spaced at round values.
#'
#' @param x A vector, matrix, data.table or list of variables to plot. If a list is provided,
#' then the bars will correspond to the mean of each sublists and individual points and
#' standard deviation will be added.
#' @param xlim Limits for the x axis.
#' @param ylim Limits for the y axis.
#' @param names.arg A vector of names to be plotted below each bar or group of bars. 
#' @param show.bar.numbers Should bar numbers be shown? If set to TRUE, the height of the bar will be
#' added, after applying the bar.numbers.FUN function. If a vector of matching length is provided, it will
#' be plotted as is.
#' @param bar.numbers.FUN A function specifying how to process bar height before plotting. 
#' If set to NULL, no value is plotted.
#' @param bar.numbers.cex Bar number labels cex.
#' @param bar.numbers.offset Bar number labels offset. Default= .5.
#' @param col Color used for bars.
#' @param beside If set to TRUE and x is a data.table or a matrix, columns are shown
#' as juxtaposed bars. Default= FALSE.
#' @param width Bar width should be smaller than 1. If beside= TRUE, then each group
#' of bar will sum to this value. Default= 0.8.
#' @param space Space between each bar. By default, it arranges them so tha they are regularly
#' spaced and aligned on increasing integer indices.
#' @param pch When x is a list, pch for the individual points.
#' @param pch.cex When x is a list, cex for the individual points.
#' @param pch.col When x is a list, color for individual points.
#' @param pch.jitter When x is a list, jitter value for individual points.
#' @param xaxt Default= "s".
#' @param yaxt Default= "s".
#' @param sd.arrow.lwd When x is a list, line width for the sd arrows. 
#' @param sd.arrow.length When x is a list, length of the sd arrow heads.
#' @param horiz Should the plot be horizontal? Default= F.
#' @param add Should the plot be added to an existing plot? Default= F.
#' @param ... Additional arguments to be passed to barplot().
#'
#' @return Plots a nice barplot with regularly spaced bars align on integers.
#'
#' @examples
#' vl_par(mfrow= c(2,2))
#' vl_barplot(c("AAAAAA"=1, "B"=2, "C"=3), horiz = F)
#' vl_barplot(c("A"=1, "B"=2, "C"=3), horiz = T)
#' vl_barplot(matrix(1:6, nrow= 2), beside = T)
#' vl_barplot(list(c(1,2,3), c(4,5,6), c(7,8,9)), beside = T, horiz= F)
#' 
#' @export
vl_barplot <- function(
    x,
    xlim= NULL,
    ylim= NULL,
    names.arg= NULL,
    ylab= NA,
    show.bar.numbers= TRUE,
    bar.numbers.FUN= function(x) round(x, .1),
    bar.numbers.cex= .7,
    bar.numbers.offset= .5,
    col= NULL,
    beside= FALSE,
    width= .8,
    space= NULL,
    pch= 16,
    pch.col= adjustcolor("grey20", .7),
    pch.cex= .5,
    pch.jitter= mean(width)/10,
    xaxt= "s",
    yaxt= "s",
    sd.arrow.lwd= .5,
    sd.arrow.length= width/10,
    horiz= F,
    add= F,
    ...
) {
  # Checks ----
  stopifnot(width<=1)
  if(is.data.table(x))
    x <- as.matrix(x)
  
  # Barplot of mean values? ----
  mean.barplot <- is.list(x)
  if(is.list(x)) {
    add.points <- x
    sd <- sapply(x, sd)
    x <- sapply(x, mean)
    if(!is.null(names(add.points)))
      names(x) <- names(add.points)
  }
  
  # Compute default names ---- 
  if(is.null(names.arg))
    names.arg <- if(is.null(colnames(x))) names(x) else colnames(x)
  
  # Compute default xlim ---- 
  if(is.null(xlim)) {
    xlim <- c(0.5, ifelse(is.null(ncol(x)), length(x), ncol(x))+.5)
  }
  
  # Compute default ylim ---- 
  if(is.null(ylim)) {
    ylim <- if(mean.barplot)
      range(unlist(add.points)) else 
        range(x)
    if(all(ylim>0))
      ylim[1] <- 0
    if(all(ylim<0))
      ylim[2] <- 0
  }
  
  # Compute space ----
  if(is.null(space)) {
    space <- rep((1-width)/width, ifelse(is.null(ncol(x)), length(x), ncol(x)))
    space[1] <- space[1]+.5
    if(beside && !is.null(nrow(x)) && nrow(x)>1) {
      space <- space*nrow(x)
      space <- c(sapply(space, function(x) c(x, 0)))
    }
  }
  
  # Adjust width when beside is TRUE ----
  if(beside && !is.null(nrow(x)) && nrow(x)>1)
    width <- width/nrow(x)
  
  # Plot barplot ----
  bar <- barplot(
    x,
    width= width,
    space= space,
    xlim= if(horiz) ylim else xlim,
    ylim= if(horiz) xlim else ylim,
    beside= beside,
    horiz= horiz,
    add= add,
    xaxt= ifelse(horiz, "s", "n"),
    yaxt= ifelse(horiz, "n", "s"),
    col= col,
    ylab= ylab,
    ...
  )
  bars <- c(bar)
  
  # Add sd and points on mean bars ----
  if(mean.barplot) {
    sapply(seq_along(bars), function(i) {
      # Standard deviation
      x.pos <- bars[i]
      segments(
        ifelse(horiz, x[i]-sd[i], x.pos), 
        ifelse(horiz, x.pos, x[i]-sd[i]),
        ifelse(horiz, x[i]+sd[i], x.pos),
        ifelse(horiz, x.pos, x[i]+sd[i]),
        lwd= sd.arrow.lwd,
        xpd= NA
      )
      segments(
        ifelse(horiz, x[i]-sd[i], x.pos-sd.arrow.length), 
        ifelse(horiz, x.pos-sd.arrow.length, x[i]-sd[i]),
        ifelse(horiz, x[i]-sd[i], x.pos+sd.arrow.length),
        ifelse(horiz, x.pos+sd.arrow.length, x[i]-sd[i]),
        lwd= sd.arrow.lwd,,
        xpd= NA
      )
      segments(
        ifelse(horiz, x[i]+sd[i], x.pos-sd.arrow.length),
        ifelse(horiz, x.pos-sd.arrow.length, x[i]+sd[i]),
        ifelse(horiz, x[i]+sd[i], x.pos+sd.arrow.length),
        ifelse(horiz, x.pos+sd.arrow.length, x[i]+sd[i]),
        lwd= sd.arrow.lwd,,
        xpd= NA
      )
      # Points
      y.pos <- add.points[[i]]
      x.pos <- jitter(rep(x.pos, length(y.pos)), amount = pch.jitter)
      if(horiz) {
        h.y.pos <- x.pos
        x.pos <- y.pos
        y.pos <- h.y.pos
      }
      points(
        x.pos,
        y.pos,
        pch= pch,
        col= pch.col,
        xpd= NA,
        cex= pch.cex
      )
    }
    )
  } else {
    if(!isFALSE(show.bar.numbers)) {
      # Compute bars height
      bar.height <- if(!beside && !is.null(nrow(x)) && nrow(x)>1)
        apply(x, sum) else
          unlist(x)
      
      # Compute labels
      bar.labels <- if(isTRUE(show.bar.numbers))
        bar.numbers.FUN(bar.height) else if (is.vector(show.bar.numbers) && length(show.bar.numbers)==length(bars))
          show.bar.numbers else
            stop("Error while computing labels to plot")
      
      # Plot 
      x.lab <- if(horiz) bar.height else bars
      y.lab <- if(horiz) bars else bar.height
      pos.lab <- if(horiz) ifelse(bar.height>0, 4, 2) else ifelse(bar.height>0, 3, 1)
      text(
        x = x.lab,
        y = y.lab,
        labels = bar.labels,
        cex = bar.numbers.cex,
        xpd= NA,
        pos= pos.lab,
        offset= bar.numbers.offset
      )
    }
  }
  
  # Add x labels ----
  if(!is.null(names.arg) && yaxt != "n") {
    if(horiz) {
      axis(
        2,
        at= seq(ifelse(is.null(ncol(x)), length(x), ncol(x))),
        labels = names.arg,
        lwd= 0
      )
    } else {
      if(max(strwidth(names.arg, cex = par("cex.axis")))>.9) {
        tiltAxis(
          x= seq(ifelse(is.null(ncol(x)), length(x), ncol(x))),
          labels = names.arg
        )
      } else {
        axis(
          1,
          at= seq(ifelse(is.null(ncol(x)), length(x), ncol(x))),
          labels = names.arg,
          lwd= 0,
          padj= -1.25
        )
      }
    }
  }
  
  # Return bars
  invisible(bar)
}