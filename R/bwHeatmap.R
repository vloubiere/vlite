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
#' @param col.palette Vector of colors for the tracks (default: rainbow palette)
#' @param xlab X axis label (default: Genomic distance)
#' @param ylab Y axis label (default: Mean signal)
#' @param cleanup.cache Logical. Whether to force recomputation of cached results (default: FALSE)
#'
#' @return Plots the average signal profile and invisibly returns NULL
#'
#' @examples
#' # Import 200 genes for example
#' path <- system.file("extdata", "Drosophila_transcripts_r6.36.bed", package = "vlite")
#' genes <- importBed(path)[5001:5200]
#'
#' # Example tracks
#' tracks <- c("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PH18_merge.bw",
#'             "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PHD11_merge.bw")
#'
#' # Parameters for a nice layout
#' vl_par()
#'
#' # Simple example with two tracks
#' bwHeatmap(bed = genes,
#'           tracks = tracks,
#'           center = "start",
#'           upstream = 5000,
#'           downstream = 5000,
#'           nbins = 251L)
#'
#' # Two tracks and two subgroups of regions
#' bwHeatmap(bed = genes,
#'           by = c(rep("group.a", 50), rep("group.b", 150)),
#'           tracks = tracks,
#'           center = "start",
#'           upstream = 5000,
#'           downstream = 5000,
#'           nbins = 251L)
#' abline(v= 0, lty= 2)
#'
#' @export
bwHeatmap <- function(bed,
                      by= NULL,
                      tracks,
                      names= function(x) gsub("(.*)[.].*$", "\\1", basename(tracks)),
                      center= "center",
                      center.name= ifelse(center!="region", "Center", "Start"),
                      region.end= "End",
                      upstream= 5000L,
                      downstream= 5000L,
                      max.cutoffs= NULL,
                      max.fun= function(x) quantile(x, .999),
                      ordering.track.idx= 1,
                      ordering.fun= function(x) mean(x),
                      nbins= if(center=="region") c(50L, 150L, 50L) else 251L,
                      ignore.strand= FALSE,
                      col.palette= c("blue", "yellow"),
                      xlab= "Genomic distance",
                      ylab= NULL,
                      cleanup.cache= FALSE)
{
  # Checks
  if(is.function(names))
    names <- names(tracks)
  if(!is.null(max.cutoffs) && length(max.cutoffs)!=length(bw))
    stop("max.cutoffs should match the number of bw files.")

  # Compute signal if necessary (cache file) ----
  signal <- bwBinnedCoverage(bed = bed,
                             tracks = tracks,
                             center = center,
                             upstream = upstream,
                             downstream = downstream,
                             nbins = nbins,
                             ignore.strand= ignore.strand)

  # Compute max values ----
  if(is.null(max.cutoffs))
    max.cutoffs <- sapply(signal, function(x) max.fun(unlist(x)))

  # Compute ordering statistics ----
  ordering.var <- if(!is.null(ordering.track.idx)) {
    apply(signal[[ordering.track.idx]], 1, ordering.fun)
  }else
    seq(nrow(bed))

  # Compute order based on groups/selected track ----
  ord <- if(is.null(by)){
    order(ordering.var, decreasing= TRUE)
  } else {
    order(by, -ordering.var)
  }

  # Plot ----
  for(i in seq(signal)) {
    # Heatmap
    vl_heatmap(signal[[i]][ord,],
               cluster.rows = if(is.null(by)) FALSE else sort(by),
               cluster.cols = FALSE,
               show.rownames = FALSE,
               zlim = c(0, max.cutoffs[i]),
               show.row.clusters = "left",
               show.colnames = FALSE,
               legend.title = "Signal")
    title(main= names[i],
          xlab= xlab,
          ylab= ylab)

    # Add x axis
    if(center!="region") {
      if(upstream>0 & downstream>0) {
        axis(1,
             at = c(1, upstream/(upstream-(-downstream))*(nbins+1), nbins),
             labels = c(-upstream, center.name, downstream))
      } else {
        axis(1,
             at = c(1, nbins),
             labels = c(-upstream, downstream))
      }
    } else {
      axis(1,
           at = c(-upstream, 0, diff(xlim)*nbins[2]/sum(nbins), xlim[2]),
           labels = c(-upstream, center.name, "End", downstream))
    }
  }
}


