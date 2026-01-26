#' Calculate Feature Coverage Across Genomic Intervals
#'
#' @description
#' A wrapper around `GenomicRanges::countOverlaps()` that computes, for each genomic range in a,
#' the number of overlapping regions in b.
#'
#' @param a Query regions in any format compatible with `importBed()`.
#' @param b Target regions in any format compatible with `importBed()`.
#' @param maxgap A single integer specifying the maximum gap allowed between 2 ranges for them to
#' be considered as overlapping. e.g., with maxgap= 0, touching regions will be considered as overlapping.
#' Default= -1L (>= 1 overlapping base).
#' @param minoverlap A single integer specifying the minimum overlap between 2 ranges for them to
#' be considered as overlapping. Default= 1L.
#' @param ignore.strand If set to FALSE, only features that are on the same strand will be counted.
#' If set to TRUE (default), overlapping features are counted irrespective of their strand.
#'
#' @return A numeric vector of length nrow(a) corresponding, for each region in a, to the number
#' of overlapping regions in b. 0 means no overlaps.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-300:+", "chr1:400-500:-", "chr2:100-200:+"))
#' b <- importBed(c("chr1:50-150:+", "chr1:200-250:+", "chr1:400-450:-", "chr1:400-450:+", "chr3:100-200:+"))
#'
#' # Count overlapping features
#' covBed(a, b)
#' covBed(a, b, ignore.strand = FALSE)
#'
#' @export
covBed <- function(a,
                   b,
                   maxgap= -1L,
                   minoverlap= 1L,
                   ignore.strand= TRUE)
{
  # Checks ----
  if(maxgap > (-1L) && minoverlap > 0L) {
    warning("Maxgap >= 0L -> minoverlap automatically set to 0L.")
    minoverlap <- 0L
  }
  
  # Import for incapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Compute coverage ----
  cov <- suppressWarnings(
    GenomicRanges::countOverlaps(
      query = GRanges(a),
      subject = GRanges(b),
      maxgap = maxgap,
      minoverlap = minoverlap,
      ignore.strand= ignore.strand
    )
  )

  # Return ----
  return(cov)
}
