#' Subset overlapping or non-overlapping Regions
#'
#' @description
#' A wrapper around ?GenomicRanges::countOverlaps that subsets the genomic ranges in a that do
#' (or do not) overlap region(s) in b.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param maxgap A single integer specifying the maximum gap allowed between 2 ranges for them to
#' be considered as overlapping. Default= -1.
#' @param minoverlap A single integer specifying the minimum overlap between 2 ranges for them to
#' be considered as overlapping. Default= 0.
#' @param invert If set to TRUE, returns regions in a that have no overlap(s) in b.
#'   If set to FALSE (default), overlapping regions are returned.
#' @param ignore.strand If set to FALSE, only overlapping features that are on the same strand are reported.
#' If set to TRUE (default), overlaps will be computed irrespective of the strand.
#'
#' @return A data.table containing:
#' \itemize{
#'   \item When invert = FALSE: regions from a that overlap region(s) in b.
#'   \item When invert = TRUE: regions from a that DO NOT overlap any region(s) in b.
#' }
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-200:+", "chr1:300-400:-", "chr1:700-800:-", "chr2:500-600:+"))
#' b <- importBed(c("chr1:100-200:+", "chr1:350-450:+", "chr1:700-800:+"))
#'
#' # Find all overlapping regions
#' intersectBed(a, b)
#' intersectBed(a, b, minoverlap = 75L)
#' intersectBed(a, b, ignore.strand= FALSE)
#'
#' # Find non-overlapping regions
#' intersectBed(a, b, invert = TRUE)
#' intersectBed(a, b, invert = TRUE, ignore.strand= FALSE)
#'
#' @export
intersectBed <- function(a,
                         b,
                         maxgap= -1L,
                         minoverlap= 0L,
                         invert= FALSE,
                         ignore.strand= TRUE)
{
  # Import bed for incapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Overlaps ----
  tmp <- rev(make.unique(c(names(a), "ov")))[1]
  assign(
    tmp,
    GenomicRanges::countOverlaps(
      query = GRanges(a),
      subject = GRanges(b),
      maxgap = maxgap,
      minoverlap = minoverlap,
      ignore.strand= ignore.strand
    )
  )

  # Non-intersecting indexes ----
  res <- if(invert) {
    a[get(tmp)==0]
  } else {
    a[get(tmp)>0]
  }

  # Return, preserving original order ----
  return(res)
}
