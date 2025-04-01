#' Find Overlapping or Non-overlapping Genomic Regions
#'
#' @description
#' Identifies genomic regions that overlap (or don't overlap) between two sets,
#' with support for strand-specific matching and minimum overlap requirements.
#' Similar to `bedtools intersect`.
#'
#' @param a Query regions in any format compatible with `importBed()`:
#' \itemize{
#'   \item Character vector of ranges ("chr:start-end[:strand]")
#'   \item GRanges object
#'   \item data.frame/data.table with required columns
#'   \item Path to a BED file
#' }
#' @param b Target regions in same format as `a`. These regions are checked
#'   for overlaps with regions in `a`.
#' @param ignore.strand Logical. If `TRUE`, finds overlaps regardless of strand.
#'   If `FALSE`, only considers overlaps between matching strands. Default is `TRUE`.
#' @param min.overlap Integer or integer vector. Minimum required overlap size in
#'   base pairs. Can be:
#' \itemize{
#'   \item Single value: Applied to all regions
#'   \item Vector: Different threshold for each region in `a`
#' }
#' Default is `1L` (any overlap).
#' @param invert Logical. If `TRUE`, returns non-overlapping regions from `a`.
#'   If `FALSE`, returns overlapping regions. Default is `FALSE`.
#'
#' @details
#' **Overlap Detection**:
#'
#' Regions overlap if they:
#' - Share the same chromosome
#' - Have overlapping coordinates ≥ `min.overlap` bases
#' - Match strands (if `ignore.strand = FALSE`)
#'
#' **Strand Handling**:
#'
#' - When `ignore.strand = TRUE`:
#'   * Finds overlaps regardless of strand
#'   * Useful for strand-independent features
#' - When `ignore.strand = FALSE`:
#'   * Only considers overlaps between matching strands
#'   * Required for strand-specific analyses
#'
#' **Minimum Overlap**:
#'
#' - Single value: Same threshold for all regions
#' - Vector: Different thresholds per region
#' - Regions with overlap < threshold are treated as non-overlapping
#'
#' @return A data.table containing:
#' \itemize{
#'   \item When `invert = FALSE`: Regions from `a` that overlap with `b`
#'   \item When `invert = TRUE`: Regions from `a` that don't overlap with `b`
#' }
#' Original column order and data types are preserved.
#'
#' @examples
#' # Create example regions
#' query <- importBed(c(
#'   "chr1:100-200:+",  # Complete overlap
#'   "chr1:300-400:-",  # Partial overlap
#'   "chr1:500-600:+",  # No overlap
#'   "chr1:700-800:-"   # Strand-specific overlap
#' ))
#'
#' target <- importBed(c(
#'   "chr1:100-200:+",  # Matches first region exactly
#'   "chr1:350-450:+",  # Overlaps second region (different strand)
#'   "chr1:700-800:+"   # Overlaps fourth region (different strand)
#' ))
#'
#' # Find all overlapping regions
#' overlaps <- intersectBed(query, target)
#' print("Overlapping regions:")
#' print(overlaps)
#'
#' # Find strand-specific overlaps
#' strand_overlaps <- intersectBed(query, target, ignore.strand = FALSE)
#' print("Strand-specific overlaps:")
#' print(strand_overlaps)
#'
#' # Find non-overlapping regions
#' non_overlaps <- intersectBed(query, target, invert = TRUE)
#' print("Non-overlapping regions:")
#' print(non_overlaps)
#'
#' # Require minimum 50bp overlap
#' min_50bp <- intersectBed(query, target, min.overlap = 50L)
#' print("Regions with ≥50bp overlap:")
#' print(min_50bp)
#'
#' # Different overlap requirements per region
#' varying_overlap <- intersectBed(query, target,
#'                                min.overlap = c(10L, 50L, 20L, 30L))
#'
#' @export
intersectBed <- function(a,
                         b,
                         ignore.strand= TRUE,
                         min.overlap= 1L,
                         invert= FALSE)
{
  # Overlap a and b ----
  idx <- overlapBed(a,
                    b,
                    ignore.strand = ignore.strand,
                    all.a = TRUE)

  # Intersection indexes ----
  inter <- unique(idx[overlap.width >= min.overlap, idx.a])

  # Non-intersecting indexes ----
  if(invert)
    inter <- setdiff(seq(nrow(a)), inter)

  # Return, preserving original order ----
  return(a[inter])
}
