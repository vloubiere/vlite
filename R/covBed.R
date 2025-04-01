#' Calculate Feature Coverage Across Genomic Intervals
#'
#' @description
#' Counts overlapping features from a target set for each interval in a query set.
#' Supports both strand-aware and strand-agnostic counting, similar to
#' `bedtools coverage`.
#'
#' @param a Query regions in any format compatible with `importBed()`:
#' \itemize{
#'   \item Character vector of ranges ("chr:start-end[:strand]")
#'   \item GRanges object
#'   \item data.frame/data.table with required columns
#'   \item Path to a BED file
#' }
#' @param b Target regions in same format as `a`. These features are counted
#'   when they overlap regions in `a`.
#' @param ignore.strand Logical. If `TRUE`, counts overlaps regardless of strand.
#'   If `FALSE`, only counts features on matching strands. Default is `TRUE`.
#'
#' @details
#' **Coverage Calculation**:
#'
#' For each region in `a`:
#' - Identifies all overlapping features from `b`
#' - Counts number of distinct overlapping features
#' - Reports 0 for regions with no overlaps
#'
#' **Strand Handling**:
#'
#' - When `ignore.strand = TRUE`:
#'   * Counts all overlapping features
#'   * Useful for strand-independent analyses
#' - When `ignore.strand = FALSE`:
#'   * Only counts features on matching strands
#'   * Required for strand-specific analyses
#'
#' **Overlap Definition**:
#'
#' Features are counted if they:
#' - Share the same chromosome
#' - Have any overlap in coordinates
#' - Match strands (if `ignore.strand = FALSE`)
#'
#' @return A numeric vector where:
#' \itemize{
#'   \item Length equals number of regions in `a`
#'   \item Each value represents number of overlapping features from `b`
#'   \item Values ≥ 0 (0 indicates no overlapping features)
#' }
#'
#' @examples
#' # Create example regions
#' query <- importBed(c(
#'   "chr1:100-300:+",  # Region with multiple overlaps
#'   "chr1:400-500:-",  # Region with strand-specific overlap
#'   "chr2:100-200:+"   # Region with no overlaps
#' ))
#'
#' features <- importBed(c(
#'   "chr1:50-150:+",    # Overlaps first region
#'   "chr1:200-250:+",   # Overlaps first region
#'   "chr1:400-450:-",   # Overlaps second region (strand-specific)
#'   "chr1:400-450:+",   # Overlaps second region (opposite strand)
#'   "chr3:100-200:+"    # No overlaps with query
#' ))
#'
#' # Count all overlapping features
#' strand_agnostic <- covBed(query, features)
#' print("Strand-agnostic coverage:")
#' print(strand_agnostic)
#'
#' # Count only same-strand features
#' strand_specific <- covBed(query, features, ignore.strand = FALSE)
#' print("Strand-specific coverage:")
#' print(strand_specific)
#'
#' # Complex example with nested features
#' complex_query <- importBed(c(
#'   "chr1:100-500:+"    # Large region containing multiple features
#' ))
#'
#' complex_features <- importBed(c(
#'   "chr1:150-200:+",   # Contained within query
#'   "chr1:300-350:+",   # Contained within query
#'   "chr1:400-600:+"    # Partially overlaps query
#' ))
#'
#' complex_coverage <- covBed(complex_query, complex_features)
#' print("Complex coverage example:")
#' print(complex_coverage)
#'
#' @export
covBed <- function(a,
                   b,
                   ignore.strand= TRUE)
{
  # Compute overlaps ----
  idx <- overlapBed(a,
                    b,
                    ignore.strand = ignore.strand,
                    all.a= TRUE)

  # Compute coverage ----
  cov <- idx[, sum(overlap.width>0), keyby= idx.a]$V1

  # Return ----
  return(cov)
}
