#' Find Overlapping or Non-overlapping Genomic Regions
#'
#' @description
#' Identifies genomic ranges in a that do overlap (or not) with any genomic ranges in b.
#'
#' @param a Query regions in any format compatible with ?importBed().
#' @param b Target regions in any format compatible with ?importBed().
#' @param min.overlap Integer specifying the minimum overlap width required. Can be a single value (that will
#' be applied to all regions) or a vector of length nrow(a). Default= 1L.
#' @param invert If set to TRUE, returns regions in a that have no overlap(s) in b.
#'   If set to FALSE (default), overlapping regions are returned.
#' @param ignore.strand If set to FALSE and strand column is provided, only overlapping features that are on
#' the same strand are reported. If set to TRUE (default), overlapping features on both strands are reported.
#'
#' @return A gr data.table containing:
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
#' intersectBed(a, b, min.overlap = 75L)
#' intersectBed(a, b, ignore.strand= FALSE)
#'
#' # Find non-overlapping regions
#' intersectBed(a, b, invert = TRUE)
#' intersectBed(a, b, invert = TRUE, ignore.strand= FALSE)
#'
#' @export
intersectBed <- function(a,
                         b,
                         min.overlap= 1L,
                         invert= FALSE,
                         ignore.strand= TRUE)
{
  # Import bed files ----
  a <- importBed(a)
  b <- importBed(b)

  # Overlap a and b and return intersecting indices ----
  inter <- if(min.overlap==1L) {
    # No need to compute overlap.width -> use covBed
    which(covBed(a, b, ignore.strand = ignore.strand)>0)
  } else {
    # Compute overlap.width using overlapBed
    idx <- overlapBed(a,
               b,
               ignore.strand = ignore.strand,
               all.a = TRUE)
    # Intersection indices
    unique(idx[overlap.width >= min.overlap, idx.a])
  }

  # Non-intersecting indexes ----
  if(invert)
    inter <- -inter

  # Return, preserving original order ----
  return(a[inter])
}
