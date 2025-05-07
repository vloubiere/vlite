#' Bin genomic regions and compute BigWig signal
#'
#' This function first bins genomic regions from a BED file, then quantifies signal from one or more BigWig files.
#' For each BigWig file, it returns a matrix where rows represent regions from the BED file and columns represent bins.
#'
#' @param bed A BED-formatted object compatible with `rtracklayer::import()`. Can be a file path or a GRanges object.
#' @param tracks Character vector of paths to BigWig (.bw) files.
#' @param center Character string specifying the anchor point for region resizing:
#'   * `center`: Resize around the region center (default)
#'   * `start`: Resize from region start
#'   * `end`: Resize from region end
#'   * `region`: Extend or contract entire region boundaries
#' @param upstream Integer specifying bases to extend upstream (default: 5000L).
#' @param downstream Integer specifying bases to extend downstream (default: 5000L).
#' @param nbins Integer vector specifying number of bins:
#'   * For `center = 'center'/'start'/'end'"'`: Length 1 vector (default: 251L)
#'   * For `center = 'region'"'`: Length 3 vector for upstream/body/downstream bins (default: c(50L, 150L, 50L))
#' @param ignore.strand Logical indicating whether to ignore strand information (default: FALSE).
#' @param cleanup.cache Logical. Whether to force recomputation of cached results (default: FALSE).
#'
#' @return A list of data.frames, one per BigWig file. Each matrix has:
#'   * Rows: Regions from the BED file
#'   * Columns: Signal values for each bin
#'   * Rownames: Region names (if present in BED)
#'   * Colnames: Bin center coordinates (and pseudocoordinates when center= 'region')
#'
#' @examples
#'
#' # Example regions (100 genes)
#' path <- system.file("extdata", "Drosophila_transcripts_r6.36.bed", package = "dtBedTools")
#' bed <- importBed(path)[5000:5100]
#'
#' # Example tracks
#' tracks <- c("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PH18_merge.bw",
#'             "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PHD11_merge.bw")
#'
#' @export
bwBinnedCoverage <- function(bed,
                             tracks,
                             center= "center",
                             upstream= 5000L,
                             downstream= 5000L,
                             nbins= if(center=="region") c(50L, 150L, 50L) else 251L,
                             ignore.strand= FALSE,
                             cleanup.cache= FALSE)
{
  # Checks ----
  if(length(upstream)>1)
    stop("'upstream' should be a unique integer value.")
  if(length(downstream)>1)
    stop("'downstream' should be a unique integer value.")

  # Use a temp directory for caching signal files ----
  cache_dir <- tempdir()
  signal.params <- list(bed,
                        tracks,
                        center,
                        upstream,
                        downstream,
                        nbins,
                        ignore.strand)
  signal.key <- digest::digest(signal.params)
  signal.cache <- file.path(cache_dir, paste0(signal.key, ".rds"))

  # If file does not exist ----
  if(cleanup.cache | !file.exists(signal.cache)) {

    # Hard copy for incapsulation ----
    regions <- importBed(bed)

    # Resize and bin ----
    bins <- if(center!="region")
    {
      # Check
      if(length(nbins)!=1L)
        stop("When center is set to 'center', 'start' or 'end', nbins should be of length 1")

      # Resize regions
      resized <- resizeBed(regions,
                           center = center,
                           upstream = upstream,
                           downstream = downstream,
                           ignore.strand = ignore.strand)

      # Compute bins
      binBed(resized,
             nbins = nbins,
             center.nbins = TRUE,
             ignore.strand = ignore.strand)

    }else if(center=="region")
    {
      if(length(nbins)!=3L)
        stop("When center is set to 'region', nbins should be of length 3")
      # Upstream
      up <- resizeBed(regions,
                      center = "start",
                      upstream = upstream,
                      downstream = 0,
                      ignore.strand = ignore.strand)
      upBins <- binBed(up,
                       nbins= nbins[1],
                       ignore.strand = ignore.strand,
                       center.nbins = TRUE)
      # Downstream
      down <- resizeBed(regions,
                        center = "end",
                        upstream = 0,
                        downstream = downstream,
                        ignore.strand = ignore.strand)
      downBins <- binBed(down,
                         nbins= nbins[3],
                         ignore.strand = ignore.strand,
                         center.nbins = TRUE)
      # Adjust middle regions after centering up/down bins'
      middle <- importBed(regions)
      middle[, start:= ifelse(up$end<down$end, up$end, down$end)+1L]
      middle[, end:= ifelse(up$start>down$start, up$start, down$start)-1L]
      # Bin middle region
      middleBins <- binBed(middle,
                           nbins= nbins[2],
                           ignore.strand = ignore.strand,
                           center.nbins = FALSE)
      # Adjust bin indexes
      middleBins[, bin.idx:= bin.idx+nbins[1]]
      downBins[, bin.idx:= bin.idx+sum(nbins[1:2])]
      # rbind
      rbind(upBins,
            middleBins,
            downBins)
    }else
      stop("center should be one of 'center', 'start', 'end' or 'region'")

    # Quantif tracks ----
    .q <- parallel::mclapply(tracks, function(x) {
      bwCoverage(bins, x)
    },
    mc.preschedule = T,
    mc.cores = max(c(1, data.table::getDTthreads()-1)))
    .q <- do.call(cbind, .q)
    .q <- as.data.table(.q)

    # Combine ----
    obj <- cbind(bins, .q)
    setnames(obj,
             make.unique(names(obj)))

    # Compute bins center ----
    bins.center <- if(center!="region") {
      seq(-upstream, downstream, length.out= nbins)
    } else {
      up <- seq(-upstream, 0, length.out= nbins[1])
      bin.width <- diff(up)[1]
      mid <- seq(0+bin.width, bin.width*nbins[2], length.out= nbins[2])
      down.start <- bin.width*(nbins[2]+1)
      down <- seq(down.start, down.start+downstream, length.out= nbins[3])
      c(up, mid, down)
    }


    # dcast ----
    cols <- data.table::last(names(obj), n = ncol(.q))
    res <- lapply(cols, function(col) {
      .c <- dcast(obj, line.idx~bin.idx, value.var = col)
      .c <- as.matrix(.c, 1)
      colnames(.c) <- as.character(bins.center)
      .c
    })

    # Save in cache ----
    saveRDS(res, signal.cache)
  } else {

    # Import signal from cache memory ----
    res <- readRDS(signal.cache)
  }

  # Return ----
  return(res)
}
