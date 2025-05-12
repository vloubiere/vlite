#' Find Overlapping Regions Between Two Sets of Genomic Intervals
#'
#' @description
#' For each genomic range in a, identifies overlaps with genomic ranges in b.
#'
#' @param a Query regions in any format compatible with ?importBed().
#' @param b Target regions in any format compatible with ?importBed().
#' @param all.a If set to FALSE, only returns the regions in a for each at least one overlap was found.
#' If set to TRUE (default), reports all regions in a, including those without overlaps (see return values).
#' @param ignore.strand If set to FALSE, only reports overlaps between regions that are on the same strand.
#' If set to TRUE (default), reports overlaps on both strands.
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item idx.a: line index of the region in a.
#'   \item idx.b: line index of the overlapping region in b (NA if no overlap).
#'   \item overlap.start: overlap starting position (NA if no overlap).
#'   \item overlap.end: overlap ending position (NA if no overlap).
#'   \item overlap.width: overlap width (0 if no overlap).
#' }
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-200:+", "chr1:300-400:-", "chr2:100-200:+"))
#' b <- importBed(c("chr1:150-250:+", "chr1:350-450:+", "chr1:500-600:+"))
#'
#' # Find all overlaps (strand-agnostic)
#' overlapBed(a, b)
#' overlapBed(a, b, all.a = FALSE)
#'
#' # Find strand-specific overlaps
#' overlapBed(a, b, ignore.strand = FALSE)
#'
#' @export
overlapBed_2 <- function(a,
                       b,
                       all.a= TRUE,
                       ignore.strand= TRUE)
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
  idx <- foverlaps(b,
                   a,
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
