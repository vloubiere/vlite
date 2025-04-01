#' Find Overlapping Regions Between Two Sets of Genomic Intervals
#'
#' @description
#' Identifies and characterizes overlaps between two sets of genomic intervals.
#' Reports overlap coordinates, widths, and can operate in both strand-aware
#' and strand-agnostic modes.
#'
#' @param a Query regions in any format compatible with `importBed()`:
#' \itemize{
#'   \item Character vector of ranges ("chr:start-end[:strand]")
#'   \item GRanges object
#'   \item data.frame/data.table with required columns
#'   \item Path to a BED file
#' }
#' @param b Target regions in same format as `a`. These are checked for
#'   overlaps with regions in `a`.
#' @param ignore.strand Logical. If `TRUE`, finds overlaps regardless of strand.
#'   If `FALSE`, only reports overlaps between regions on matching strands.
#'   Default is `TRUE`.
#' @param all.a Logical. Controls output comprehensiveness:
#' \itemize{
#'   \item `TRUE`: Reports all regions from `a`, including those without overlaps (default)
#'   \item `FALSE`: Reports only regions from `a` that have overlaps in `b`
#' }
#'
#' @details
#' **Overlap Detection**:
#'
#' Regions overlap if they:
#' - Share the same chromosome
#' - Have overlapping coordinates
#' - Match strands (if `ignore.strand = FALSE`)
#'
#' **Overlap Calculations**:
#'
#' For overlapping regions:
#' - Start = max(region_a_start, region_b_start)
#' - End = min(region_a_end, region_b_end)
#' - Width = end - start + 1
#'
#' **Strand Handling**:
#'
#' - When `ignore.strand = TRUE`:
#'   * Overlaps found regardless of strand
#'   * Useful for strand-independent features
#' - When `ignore.strand = FALSE`:
#'   * Only reports overlaps on matching strands
#'   * Required for strand-specific analyses
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item idx.a: Integer index of region in set A
#'   \item idx.b: Integer index of overlapping region in set B (NA if no overlap)
#'   \item overlap.start: Start position of overlap (NA if no overlap)
#'   \item overlap.end: End position of overlap (NA if no overlap)
#'   \item overlap.width: Width of overlap (0 if no overlap)
#' }
#'
#' @examples
#' # Create example regions
#' a <- importBed(c(
#'   "chr1:100-200:+",  # Region A1
#'   "chr1:300-400:-",  # Region A2
#'   "chr2:100-200:+"   # Region A3
#' ))
#'
#' b <- importBed(c(
#'   "chr1:150-250:+",  # Overlaps A1
#'   "chr1:350-450:-",  # Overlaps A2
#'   "chr1:500-600:+"   # No overlaps
#' ))
#'
#' # Find all overlaps (strand-agnostic)
#' overlaps <- overlapBed(a, b)
#'
#' # Show overlaps with source regions
#' cbind(overlaps,
#'       "a_region" = a[overlaps$idx.a],
#'       "b_region" = b[overlaps$idx.b])
#'
#' # Find strand-specific overlaps
#' strand_overlaps <- overlapBed(a, b, ignore.strand = FALSE)
#'
#' # Only report regions with overlaps
#' matching_only <- overlapBed(a, b, all.a = FALSE)
#'
#' # Complex overlap example
#' a2 <- importBed(c(
#'   "chr1:100-300:+",  # Multiple overlaps
#'   "chr2:200-400:-"   # No overlaps
#' ))
#'
#' b2 <- importBed(c(
#'   "chr1:50-150:+",   # Partial overlap
#'   "chr1:250-350:+",  # Partial overlap
#'   "chr1:175-225:+"   # Contained within
#' ))
#'
#' complex_overlaps <- overlapBed(a2, b2)
#'
#' @export
overlapBed <- function(a,
                       b,
                       ignore.strand= TRUE,
                       all.a= TRUE)
{
  # Hard copy for incapsulation ----
  a <- importBed(a)
  b <- importBed(b)

  # Add index columns ----
  a[, idx.a:= .I]
  b[, idx.b:= .I]

  # Check ----
  .cols <- if(!ignore.strand && "strand" %in% names(a) && "strand" %in% names(b))
    c("seqnames", "strand", "start", "end") else
      c("seqnames", "start", "end")

  # Key bed files ----
  setkeyv(a, .cols)
  setkeyv(b, .cols)

  # foverlaps ----
  idx <- foverlaps(a,
                   b,
                   nomatch= if(all.a) NA else NULL)

  # Compute overlap.start, overlap.end, overlap.width
  idx[, overlap.start:= ifelse(i.start>start, i.start, start)]
  idx[, overlap.end:= ifelse(i.end<end, i.end, end)]
  idx[, overlap.width:= overlap.end-overlap.start+1]
  idx[is.na(overlap.width), overlap.width:= 0]

  # Sort rows and select columns ----
  setorderv(idx, c("idx.a", "idx.b"))
  idx <- idx[, .(idx.a, idx.b, overlap.start, overlap.end, overlap.width)]

  # Return ----
  return(idx)
}
