#' Subtract Genomic Regions from Reference Intervals
#'
#' @description
#' Creates new genomic intervals by removing overlapping regions from reference
#' intervals, similar to `bedtools subtract`. Supports strand-specific operations
#' and minimum width filtering.
#'
#' @param a Reference regions in any format compatible with `importBed()`:
#' \itemize{
#'   \item Character vector of ranges ("chr:start-end[:strand]")
#'   \item GRanges object
#'   \item data.frame/data.table with required columns
#'   \item Path to a BED file
#' }
#' @param b Regions to subtract from `a`, in same format as `a`.
#' @param min.width Integer. Minimum width required for resulting regions.
#'   Regions smaller than this are discarded. Default is `1L`.
#' @param ignore.strand Logical. If `TRUE`, subtracts overlaps regardless of strand.
#'   If `FALSE`, only subtracts regions on matching strands. Default is `TRUE`.
#'
#' @details
#' **Subtraction Process**:
#'
#' 1. Overlapping regions in `b` are first merged
#' 2. For each region in `a`:
#'    - Identifies overlaps with merged `b` regions
#'    - Removes overlapping portions
#'    - Creates new regions from remaining segments
#'    - Filters by minimum width
#'
#' **Strand Handling**:
#'
#' - When `ignore.strand = TRUE`:
#'   * Subtracts all overlapping regions
#'   * Useful for strand-independent features
#' - When `ignore.strand = FALSE`:
#'   * Only subtracts regions on matching strands
#'   * Preserves strand-specific information
#'
#' **Width Filtering**:
#'
#' - Regions are filtered after subtraction
#' - Only regions ≥ `min.width` are retained
#' - Prevents creation of very small fragments
#'
#' @return A data.table containing:
#' \itemize{
#'   \item Remaining portions of regions from `a`
#'   \item All original columns from `a` preserved
#'   \item Coordinates adjusted to reflect subtraction
#'   \item Only regions meeting minimum width requirement
#' }
#'
#' @examples
#' # Create example regions
#' reference <- importBed(c(
#'   "chr1:100-500:+",    # Large region to subtract from
#'   "chr1:700-900:-",    # Region with strand
#'   "chr2:100-300:+"     # No overlaps
#' ))
#'
#' to_subtract <- importBed(c(
#'   "chr1:200-300:+",    # Internal region
#'   "chr1:400-600:+",    # Overlaps end
#'   "chr1:700-800:-",    # Same strand as second region
#'   "chr1:700-800:+"     # Different strand from second region
#' ))
#'
#' # Basic subtraction (strand-agnostic)
#' result <- subtractBed(reference, to_subtract)
#' print("Regions after subtraction:")
#' print(result)
#'
#' # Strand-specific subtraction
#' strand_result <- subtractBed(reference, to_subtract,
#'                             ignore.strand = FALSE)
#' print("Strand-specific subtraction:")
#' print(strand_result)
#'
#' # Require minimum width of 50bp
#' min_width_result <- subtractBed(reference, to_subtract,
#'                                min.width = 50L)
#' print("Regions ≥50bp after subtraction:")
#' print(min_width_result)
#'
#' # Complex example with multiple overlaps
#' complex_ref <- importBed("chr1:100-1000:+")
#' complex_sub <- importBed(c(
#'   "chr1:200-300:+",    # Creates gap in middle
#'   "chr1:400-500:+",    # Creates second gap
#'   "chr1:900-1100:+"    # Overlaps end
#' ))
#'
#' complex_result <- subtractBed(complex_ref, complex_sub)
#' print("Complex subtraction result:")
#' print(complex_result)
#'
#' @export
subtractBed <- function(a,
                        b,
                        min.width= 1L,
                        ignore.strand= TRUE)
{
  # Hard copy for incapsulation ----
  a <- importBed(a)

  # Collapse features to subtract ----
  coll <- collapseBed(b,
                      ignore.strand = ignore.strand)

  # Compute overlaps ----
  ov <- overlapBed(a,
                   coll,
                   ignore.strand = ignore.strand,
                   all.a = TRUE)
  # Subtract ----
  sub <- ov[, {
    .(seqnames= a$seqnames[idx.a],
      start= na.omit(c(a$start[idx.a], overlap.end+1)),
      end= na.omit(c(overlap.start-1, a$end[idx.a])))
  }, idx.a]

  # Select subtracted regions larger than min.width ----
  sub <- sub[end-start+1 >= min.width]

  # cbind `a` and remove index ----
  sub <- cbind(sub[, !"idx.a"],
               a[sub$idx.a, !c("seqnames", "start", "end")])

  # Return ----
  return(sub)
}
