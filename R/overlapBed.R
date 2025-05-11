#' Find Overlapping Regions Between Two Sets of Genomic Intervals
#'
#' @description
#' For each genomic range in a, identifies overlaps with genomic ranges in b.
#'
#' @param a Query regions in any format compatible with ?importBed().
#' @param b Target regions in any format compatible with ?importBed().
#' @param all.a If set to FALSE, only returns the regions in a for each at least one overlap was found.
#' If set to TRUE (default), reports all regions in a, including those without overlaps (see return values).
#' @param ignore.strand If set to FALSE and strand column is provided, only reports overlaps between regions
#' that are on the same strand. If set to TRUE (default), reports overlaps on both strands.
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
overlapBed <- function(a,
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

  # columns to overlap ----
  .cols <- if(!ignore.strand && "strand" %in% names(a) && "strand" %in% names(b))
    c("seqnames", "start<=end", "end>=start", "strand") else
      c("seqnames", "start<=end", "end>=start")

  # Compute overlaps ----
  idx.b <- b[a, .(.(idx.b)), by= .EACHI, on= .cols]$V1
  idx <- data.table(idx.a= rep(a$idx.a, lengths(idx.b)),
                               idx.b= unlist(idx.b))

  # Remove regions in a with no overlaps in b ----
  if(!all.a)
    idx <- idx[!is.na(idx.b)]

  # Compute overlap.start, overlap.end, overlap.width
  coor <- cbind(a[idx$idx.a, .(start, end)],
                b[idx$idx.b, .(b.start= start, b.end= end)])
  coor[, overlap.start:= ifelse(start>b.start, start, b.start)]
  coor[, overlap.end:= ifelse(end<b.end, end, b.end)]
  coor[, overlap.width:= overlap.end-overlap.start+1]
  coor[is.na(overlap.width), overlap.width:= 0]

  # Final result ----
  res <- cbind(idx,
               coor[, .(overlap.start, overlap.end, overlap.width)])

  # Return ----
  return(res)
}
