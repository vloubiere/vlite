#' Plot Average BigWig Signal Tracks
#'
#' Creates a visualization of average signal profiles from BigWig files around genomic features,
#' with optional grouping and customizable appearance.
#'
#' @param bed A GRanges or data.frame object containing genomic regions of interest
#' @param by Optional factor for grouping the data (e.g., by condition or feature type)
#' @param tracks Character vector of paths to BigWig files
#' @param names Function or character vector to generate track names. Default extracts names from file paths
#' @param center Method for centering profiles: "center" for point-centric or "region" for meta-gene like plots
#' @param center.name Label for the center point in the plot (default: "center" or "Start", when center=="region")
#' @param region.end When center is set to region, label for the end of the region. Default= "End".
#' @param upstream Integer. Distance upstream of center to plot (default: 5000)
#' @param downstream Integer. Distance downstream of center to plot (default: 5000)
#' @param nbins Integer vector. Number of bins for signal averaging. For center="region", expects c(upstream_bins, body_bins, downstream_bins). For center="center", single integer (default: 251)
#' @param ignore.strand Logical. If TRUE, ignores strand information (default: FALSE)
#' @param plot Should the plot be generated? Default= TRUE.
#' @param add Logical. If TRUE, adds to existing plot (default: FALSE)
#' @param col.palette Vector of colors for the tracks (default: rainbow palette)
#' @param adj.mean.color Numeric. Adjustment factor for mean line color transparency (default: 0.8)
#' @param adj.se.color Numeric. Adjustment factor for standard error area transparency (default: 0.3)
#' @param xlim Numeric vector. X-axis limits (default: NULL, automatically determined)
#' @param ylim Numeric vector. Y-axis limits (default: NULL, automatically determined)
#' @param xlab X axis label (default: Genomic distance)
#' @param ylab Y axis label (default: Mean signal)
#' @param legend Logical. Whether to display legend (default: TRUE)
#' @param legend.x x position for legend. Default= par("usr")[2].
#' @param legend.y y position for legend. Default= par("usr")[4].
#' @param legend.pos Character. Legend position (default: "topright")
#' @param legend.cex Numeric. Legend text size (default: 0.7)
#' @param cleanup.cache Logical. Whether to force recomputation of cached results (default: FALSE)
#'
#' @return Plots the average signal profile and invisibly returns NULL
#'
#' @examples
#' # Sample 200 random genes
#' all_genes <- GenomicFeatures::genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#' set.seed(42)
#' random_genes <- all_genes[sample(length(all_genes), 200)]
#'
#' # Example tracks
#' tracks <- c("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PH18_merge.bw",
#'             "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PHD11_merge.bw")
#'
#' # Parameters for a nice layout
#' vl_par()
#'
#' # Simple example with two tracks
#' meanSignal <- bwAverageTrack(bed = random_genes,
#'                              tracks = tracks,
#'                              center = "start",
#'                              upstream = 5000,
#'                              downstream = 5000,
#'                              nbins = 251L)
#'
#' # Two tracks and two subgroups of regions
#' meanSignal <- bwAverageTrack(bed = random_genes,
#'                              by = c(rep("group.a", 50), rep("group.b", 150)),
#'                              tracks = tracks,
#'                              center = "start",
#'                              upstream = 5000,
#'                              downstream = 5000,
#'                              nbins = 251L)
#' abline(v= 0, lty= 2)
#'
#' # Example with anchored region (Genomic distance is now a pseudodistance)
#' meanSignal <- bwAverageTrack(bed = random_genes,
#'                              tracks = tracks,
#'                              center = "region",
#'                              upstream = 5000,
#'                              downstream = 5000)
#' abline(v= 0, lty= 2)
#'
#' @export
bwAverageTrack <- function(bed,
                           by= NULL,
                           tracks,
                           names= function(x) gsub("(.*)[.].*$", "\\1", basename(tracks)),
                           center= "center",
                           center.name= ifelse(center!="region", "Center", "Start"),
                           region.end= "End",
                           upstream= 5000L,
                           downstream= 5000L,
                           nbins= if(center=="region") c(50L, 150L, 50L) else 251L,
                           ignore.strand= FALSE,
                           plot= TRUE,
                           add= FALSE,
                           col.palette= rainbow(7)[-7],
                           adj.mean.color= .8,
                           adj.se.color= .3,
                           xlim= NULL,
                           ylim= NULL,
                           xlab= "Genomic distance",
                           ylab= "Mean signal",
                           legend= TRUE,
                           legend.x= par("usr")[2],
                           legend.y= par("usr")[4],
                           legend.cex= .7,
                           cleanup.cache= FALSE)
{
  # Checks ----
  if(is.function(names))
    names <- names(tracks)
  if(!is.null(by))
    names <- paste0(rep(names, each= length(unique(by))), " @ ", sort(unique(by)))

  # Compute signal if necessary (cache file) ----
  signal <- bwBinnedCoverage(bed = bed,
                             tracks = tracks,
                             center = center,
                             upstream = upstream,
                             downstream = downstream,
                             nbins = nbins,
                             ignore.strand= ignore.strand,
                             cleanup.cache = cleanup.cache)
  signal <- lapply(signal, as.data.table)

  # Split groups using by ----
  if(!is.null(by)) {
    signal <- lapply(signal, function(x) split(x, by))
    signal <- unlist(signal, recursive = FALSE, use.names = FALSE)
  }

  # Compute mean signal, standard error, bin.x.pos and x.pos ----
  Nregions <- sapply(signal, nrow)
  signalMean <- lapply(signal, apply, 2, mean, na.rm= TRUE)
  signalSe <- lapply(signal, function(x) {
    apply(x, 2, function(y) {
      sd(y, na.rm = TRUE)/sqrt(length(na.omit(y)))
    })
  })
  bin.x.pos <- lapply(signal, function(x) as.numeric(colnames(x)))

  # Format data table ----
  dat <- data.table(Nregions,
                    signalMean,
                    signalSe,
                    bin.x.pos)

  # Add names and colors ----
  dat[, name:= paste0(names, " (n= ", formatC(Nregions, big.mark = ","), ")")]
  dat[, col:= colorRampPalette(col.palette)(.NGRP)[.GRP], name]

  # Melt data ----
  dat <- dat[, lapply(.SD, unlist), .(name, col)]

  # Plot ----
  setattr(dat, "class", c("bwAverageTrack", "data.table", "data.frame"))
  if(plot) {
    plot.bwAverageTrack(dat,
                        center= center,
                        upstream= upstream,
                        downstream= downstream,
                        add= add,
                        col.palette= col.palette,
                        adj.mean.color= adj.mean.color,
                        adj.se.color= adj.se.color,
                        xlim= xlim,
                        ylim= ylim,
                        xlab= xlab,
                        ylab= ylab,
                        center.name= center.name,
                        legend= legend,
                        legend.x= legend.x,
                        legend.y= legend.y,
                        legend.cex= legend.cex)
  }

  # Return mean signal ----
  invisible(dat)
}


