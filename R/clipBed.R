#' Clip Genomic Regions to Defined Limits
#'
#' @description
#' Clips genomic regions in `a` to the boundaries defined by regions in `b`.
#' This function is useful for restricting regions in `a` to valid genomic
#' boundaries, such as chromosome limits or other predefined regions.
#'
#' @param a Query regions in any format compatible with `importBed()`:
#' \itemize{
#'   \item Character vector of ranges ("chr:start-end[:strand]")
#'   \item GRanges object
#'   \item data.frame/data.table with required columns
#'   \item Path to a BED file
#' }
#' @param b Target regions in same format as `a`. These regions define the
#'   boundaries for clipping regions in `a`.
#' @param min.width Integer. Minimum width required for resulting clipped regions.
#'   Regions smaller than this are discarded. Default is `1L`.
#' @param ignore.strand Logical. If `TRUE`, clips regions regardless of strand.
#'   If `FALSE`, only clips regions on matching strands. Default is `TRUE`.
#'
#' @details
#' **Clipping Process**:
#'
#' 1. Overlapping regions in `b` are first merged.
#' 2. For each region in `a`:
#'    - Identifies overlaps with merged `b` regions.
#'    - Clips `a` regions to the boundaries of overlapping `b` regions.
#'    - Discards regions smaller than `min.width`.
#'
#' **Strand Handling**:
#'
#' - When `ignore.strand = TRUE`:
#'   * Clips regions regardless of strand.
#'   * Useful for strand-independent analyses.
#' - When `ignore.strand = FALSE`:
#'   * Only clips regions on matching strands.
#'   * Preserves strand-specific information.
#'
#' **Width Filtering**:
#'
#' - Regions are filtered after clipping.
#' - Only regions ≥ `min.width` are retained.
#' - Prevents creation of very small fragments.
#'
#' @return A data.table containing:
#' \itemize{
#'   \item Clipped portions of regions from `a`.
#'   \item All original columns from `a` preserved.
#'   \item Coordinates adjusted to reflect clipping.
#'   \item Only regions meeting minimum width requirement.
#' }
#'
#' @examples
#' # Create example regions
#' query <- importBed(c(
#'   "chr1:100-500:+",    # Large region to clip
#'   "chr1:700-900:-",    # Region with strand
#'   "chr2:100-300:+"     # No overlaps
#' ))
#'
#' boundaries <- importBed(c(
#'   "chr1:200-300:+",    # Internal region
#'   "chr1:400-600:+",    # Overlaps end
#'   "chr1:700-800:-",    # Same strand as second region
#'   "chr1:700-800:+"     # Different strand from second region
#' ))
#'
#' # Basic clipping (strand-agnostic)
#' result <- clipBed(query, boundaries)
#' print("Clipped regions:")
#' print(result)
#'
#' # Strand-specific clipping
#' strand_result <- clipBed(query, boundaries,
#'                         ignore.strand = FALSE)
#' print("Strand-specific clipping:")
#' print(strand_result)
#'
#' # Require minimum width of 50bp
#' min_width_result <- clipBed(query, boundaries,
#'                            min.width = 50L)
#' print("Regions ≥50bp after clipping:")
#' print(min_width_result)
#'
#' # Complex example with multiple overlaps
#' complex_query <- importBed("chr1:100-1000:+")
#' complex_boundaries <- importBed(c(
#'   "chr1:200-300:+",    # Creates gap in middle
#'   "chr1:400-500:+",    # Creates second gap
#'   "chr1:900-1100:+"    # Overlaps end
#' ))
#'
#' complex_result <- clipBed(complex_query, complex_boundaries)
#' print("Complex clipping result:")
#' print(complex_result)
#'
#' @export
clipBed <- function(a,
                    b,
                    min.width= 1,
                    ignore.strand= TRUE)
{
  # Hard copy for incapsulation ----
  a <- importBed(a)

  # Collapse features to subtract ----
  coll <- collapseBed(b,
                      ignore.strand = ignore.strand)

  # Compute overlaps ----
  clip <- overlapBed(a,
                     coll,
                     ignore.strand = ignore.strand,
                     all.a = FALSE)

  # Overlaps coordinates correspond to clipped coor ----
  setnames(clip,
           c("overlap.start", "overlap.end"),
           c("start", "end"))

  # Select clipped regions larger than min.width ----
  clip <- clip[end-start+1 >= min.width]

  # cbind `a` ----
  clip <- cbind(a[clip$idx.a, .(seqnames)],
                clip[, .(start, end)],
                a[clip$idx.a, !c("seqnames", "start", "end")])

  # Return ----
  return(clip)
}
