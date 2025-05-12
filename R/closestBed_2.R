#' Find Closest Genomic Features Between Two Sets of Regions
#'
#' @description
#' For each genomic range in a, identifies the closest genomic ranges in b.
#'
#' @param a Query regions in any format compatible with ?importBed().
#' @param b Target regions in any format compatible with ?importBed().
#' @param n Report all features at the n-th closest distances, including ties (e.g., if three features are
#' at distance 0 and two are at distance 1000, n = 1 returns the three at distance 0; n = 2 also returns the
#' two at distance 1000). Default = 1L.
#' @param min.dist Minimum absolute distance allowed between features. Features closer than this are excluded.
#' Default= 0L, meaning that both overlapping (dist = 0) and non-overlapping (dist > 0) features are reported.
#' @param max.dist Maximum absolute distance allowed between features. Features farther than this are excluded.
#' Default= Inf.
#' @param ignore.strand If set to FALSE and strand column is provided, only closest features that are on the
#' same strand are reported. If set to TRUE (default), closest features on both strands are reported. See details
#' about how the sign of the computed distance reflects the repesctive location of a and b features.
#'
#' @details
#' **Distance Calculation**:
#' - For overlapping features, a distance of 0 is returned.
#' - For non-overlapping features, the genomic distance between their closest boundaries are returned:
#'   * Negative distances indicate upstream features
#'   * Positive distances indicate downstream features
#' - If strand information is missing for features in a:
#'   * Negative distances indicate features in 5'
#'   * Positive distances indicate features in 3'
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item idx.a: line index of the regions in a.
#'   \item idx.b: line index of the closest region(s) in b.
#'   \item dist: distance between closest regions (see details).
#' }
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr2R:100-200:+",
#'                  "chr2R:100-200:-"))
#' b <- importBed(c("chr2R:300-400:+",
#'                  "chr2R:300-400:-",
#'                  "chr2R:3000-3100:+"))
#'
#' # Find single closest features
#' closestBed(a, b)[] # Irrespective of the strand
#' closestBed(a, b, ignore.strand = F)[] # Only consider features that are on the same strand
#' closestBed(a[, !"strand"], b)[] # No strand is provided; a features are considered as '+'
#'
#' # Return all features at the second closest distance (including ties):
#' closestBed(a, b, n= 2)[]
#'
#' @export
closestBed_2 <- function(a,
                         b,
                         n= 1L,
                         min.dist= 0L,
                         max.dist= Inf,
                         ignore.strand= TRUE)
{
  # Checks ----
  if(min.dist<0L | max.dist<0L)
    stop("min.dist and max.dist, should be positive integers")

  # Import bed files ----
  a <- GRanges(a)
  b <- GRanges(b)

  # Compute pairwise distances ----
  lapply(seq_along(a), function(i) GenomicRanges::distance(a[i], b))
}
