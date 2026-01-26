#' Clip Genomic Regions to Defined Limits
#'
#' @description
#' A wrapper around ?GenomicRanges::findOverlaps that clips the regions in a that extend beyond
#' the regions in b.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param ignore.strand If set to FALSE and strand column is provided, only clips boundaries
#' that are on the same strand. If set to TRUE (default), clips boundaries regardless of their strand.
#'
#' @return A data.table containing the remaining portions of a, after clipping them using the boundaries
#' defined in b.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-500:+"))
#' b <- importBed(c("chr1:200-300:+", "chr1:400-450:-", "chr1:425-475:-"))
#'
#' # Basic example
#' clipBed(a, b)[]
#'
#' # Strand-specific
#' clipBed(a, b, ignore.strand = FALSE)[]
#'
#' @export
clipBed <- function(a,
                    b,
                    ignore.strand= TRUE)
{
  # Import for incapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Collapse regions in b ----
  coll <- collapseBed(b, ignore.strand= ignore.strand)

  # Compute overlaps with a ----
  ov <- vlite::overlapBed(
    a = a,
    b = coll,
    ignore.strand = ignore.strand
  )

  # Resize overlaps ----
  idx.a <- ov$idx.a
  res <- a[(idx.a), env= list(idx.a= idx.a)]
  res$start <- ov$overlap.start
  res$end <- ov$overlap.end

  # Return ----
  return(res)
}
