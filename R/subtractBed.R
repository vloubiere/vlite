#' Subtract Genomic Regions from Reference Intervals
#'
#' @description
#' A wrapper around ?GenomicRanges::setdiff that subtracts the regions in b to the genomic ranges in a.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param minoverlap A single integer specifying the minimum overlap with regions in for them to
#' be subtracted. Default= 1L.
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
                        minoverlap= 1L,
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
      minoverlap= minoverlap,
      ignore.strand= ignore.strand
    )
  )

  # Retrieve row indices from a ----
  tmp <- rev(make.unique(c(names(a), "idx.a")))[1]
  assign(
    tmp,
    rep(seq(nrow(a)), lengths(sub))
  )

  # Subtracted result ----
  res <- cbind(
    as.data.table(unlist(sub))[, !"width"],
    a[get(tmp), !c("seqnames", "start", "end", "strand")]
  )

  # Return ----
  return(res)
}
