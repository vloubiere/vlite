#' Find Overlapping Regions Between Two Sets of Genomic Intervals
#'
#' @description
#' A wrapper around ?GenomicRanges::findOverlaps() that identifies, for each genomic range in a,
#' the overlapping regions in b.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param maxgap A single integer specifying the maximum gap allowed between 2 ranges for them to
#' be considered as overlapping. Default= -1.
#' @param minoverlap A single integer specifying the minimum overlap between 2 ranges for them to
#' be considered as overlapping. Default= 0.
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
#' Of note, regions in a with no overlap in b are not returned.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-200:+", "chr1:300-400:-", "chr2:100-200:+"))
#' b <- importBed(c("chr1:150-250:+", "chr1:350-450:+", "chr1:500-600:+"))
#'
#' # Find all overlaps
#' overlapBed(a, b)
#'
#' # Find strand-specific overlaps
#' overlapBed(a, b, ignore.strand = FALSE)
#'
#' @export
overlapBed <- function(a,
                       b,
                       maxgap= -1L,
                       minoverlap= 0L,
                       ignore.strand= TRUE)
{
  # Import bed for incapsulation ----
  a <- vlite::importBed(a)[, c("seqnames", "start", 'end', "strand")]
  b <- vlite::importBed(b)[, c("seqnames", "start", 'end', "strand")]

  # Overlaps ----
  ov <- suppressWarnings(
    GenomicRanges::findOverlaps(
      query = GenomicRanges::GRanges(a),
      subject = GenomicRanges::GRanges(b),
      maxgap = maxgap,
      minoverlap = minoverlap,
      ignore.strand = ignore.strand
    )
  )
  ov <- as.data.table(ov)
  setnames(ov, c("idx.a", "idx.b"))

  # Compute overlap.start, overlap.end, overlap.width
  coor <- cbind(a[ov$idx.a, .(start, end)],
                b[ov$idx.b, .(b.start= start, b.end= end)])
  coor[, overlap.start:= ifelse(start>b.start, start, b.start)]
  coor[, overlap.end:= ifelse(end<b.end, end, b.end)]
  coor[, overlap.width:= overlap.end-overlap.start+1]
  coor[is.na(overlap.width), overlap.width:= 0]

  # Final result ----
  res <- cbind(ov,
               coor[, .(overlap.start, overlap.end, overlap.width)])

  # Return ----
  return(res)
}
