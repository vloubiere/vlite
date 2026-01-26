#' Subtract Genomic Regions from Reference Intervals
#'
#' @description
#' A wrapper around ?GenomicRanges::setdiff that subtracts the regions in b to the genomic ranges in a.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param ignore.strand If set to FALSE, only subtracts features that are on the same strand.
#' If set to TRUE (default), subtracts overlapping feature(s) regardless of their strand(s).
#'
#' @return A data.table containing the remaining portions of a, after subtracting the regions
#' defined in b.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:200-300:-", "chr1:100-500:-"))
#' b <- importBed(c("chr1:200-300:+", "chr1:400-450:-", "chr1:425-475:-"))
#'
#' # Basic example
#' subtractBed(a, b)
#'
#' # Only subtract regions with similar strand
#' subtractBed(a, b, ignore.strand= FALSE)
#'
#' @export
subtractBed <- function(a,
                        b,
                        ignore.strand= TRUE)
{
  # Import for incapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Subtract features ----
  sub <- suppressWarnings(
    GenomicRanges::subtract(
      GenomicRanges::GRanges(a),
      GenomicRanges::GRanges(b),
      minoverlap= 1L,
      ignore.strand= ignore.strand
    )
  )
  
  # Store row indices ----
  row.idx <- rep(seq(nrow(a)), lengths(sub))

  # Make data.table ----
  sub <- as.data.table(unlist(sub))[, !"width"]
  
  # Retrieve additional columns ----
  add <- a[, !c("seqnames", "start", "end", "strand")]
  add <- a[(row.idx), env= list(row.idx= row.idx)]
  res <- cbind(sub, add)

  # Return ----
  return(res)
}
