#' Plot bwAverageTrack object
#'
#' Plotting method for bigwig average tracks
#'
#' @param dat A bwAverageTrack object as outputed by ?bwAverageTrack.
#' @param add Logical. If TRUE, adds to existing plot (default: FALSE)
#' @param col.palette Vector of colors for the tracks (default: rainbow palette)
#' @param adj.mean.color Numeric. Adjustment factor for mean line color transparency (default: 0.8)
#' @param adj.se.color Numeric. Adjustment factor for standard error area transparency (default: 0.3)
#' @param xlim Numeric vector. X-axis limits (default: NULL, automatically determined)
#' @param ylim Numeric vector. Y-axis limits (default: NULL, automatically determined)
#' @param xlab X axis label (default: Genomic distance)
#' @param ylab Y axis label (default: Mean signal)
#' @param center.name Label for the center point in the plot (default: "center" or "Start", when center=="region")
#' @param legend Logical. Whether to display legend (default: TRUE)
#' @param legend.x x position for legend. Default= par("usr")[2].
#' @param legend.y y position for legend. Default= par("usr")[4].
#' @param legend.pos Character. Legend position (default: "topright")
#' @param legend.cex Numeric. Legend text size (default: 0.7)
#' @param cleanup.cache Logical. Whether to force recomputation of cached results (default: FALSE)
#'
#' @return Plots the average signal profile.
#'
#' @export
plot.bwAverageTrack <- function(dat,
                                center,
                                upstream,
                                downstream,
                                add= FALSE,
                                col.palette= rainbow(7)[-7],
                                adj.mean.color= .8,
                                adj.se.color= .3,
                                xlim= NULL,
                                ylim= NULL,
                                xlab= "Genomic distance",
                                ylab= "Mean signal",
                                center.name= "Center",
                                legend= TRUE,
                                legend.x= par("usr")[2],
                                legend.y= par("usr")[4],
                                legend.cex= .7)
{
  # Initiate plot ----
  if(!add)
  {
    if(is.null(xlim))
      xlim <- range(dat$bin.x.pos)
    if(is.null(ylim))
      ylim <- range(dat[, c(signalMean-signalSe, signalMean+signalSe)])
    plot(NA,
         xlim= xlim,
         ylim= ylim,
         type= "n",
         xlab= NA,
         xaxt= "n",
         ylab= NA)
    if(!is.na(xlab))
      title(xlab = xlab, line = .65)
    if(!is.na(ylab))
      title(ylab = ylab, line = 1.5)
    # Add x axis
    if(center!="region") {
      axis(1,
           at = c(-upstream, 0, downstream),
           labels = c(-upstream, center.name, downstream),
           xpd= NA,
           padj= -1.25)
    } else {
      axis(1,
           at = unlist(bin.x.pos)[c(1, nbins[1], sum(nbins[1:2])+1, sum(nbins))],
           labels = c(-upstream, center.name, "End", downstream),
           xpd= NA,
           padj= -1.25)
    }
  }

  # Plot standard error and mean values ----
  dat[, {
    # SE
    polygon(c(bin.x.pos, rev(bin.x.pos)),
            c(signalMean-signalSe, rev(signalMean+signalSe)),
            border= NA,
            col= adjustcolor(col[1], adj.se.color))
    # Mean
    lines(bin.x.pos,
          signalMean[1:.N],
          col= adjustcolor(col[1], adj.mean.color))
  }, .(name, col)]

  # Legend ----
  if(legend)
    vl_legend(x= legend.x,
              y= legend.y,
              legend= unique(dat$name),
              text.col= unique(dat$col),
              cex = legend.cex,
              bty= "n",
              xpd= NA)
}
