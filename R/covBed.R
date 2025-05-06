#' Calculate Feature Coverage Across Genomic Intervals
#'
#' @description
#' For each genomic range in a, counts the number of overlapping genomic ranges in b.
#'
#' @param a Query regions in any format compatible with ?importBed().
#' @param b Target regions in any format compatible with ?importBed().
#' @param ignore.strand If set to FALSE, only features that are on the same strand will be counted.
#' If set to TRUE (default), overlapping feature on both strands are counted.
#'
#' @return A numeric vector of length nrow(a) corresponding, for each region in a, to the number
#' of overlapping regions in b. 0 means no overlaps.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-300:+", "chr1:400-500:-", "chr2:100-200:+"))
#' b <- importBed(c("chr1:50-150:+", "chr1:200-250:+", "chr1:400-450:-", "chr1:400-450:+", "chr3:100-200:+"))
#'
#' # Count overlapping features
#' covBed(a, b)
#' covBed(a, b, ignore.strand = FALSE)
#'
#' @export
covBed <- function(a,
                   b,
                   ignore.strand= TRUE)
{
  # Compute overlaps ----
  idx <- overlapBed(a,
                    b,
                    ignore.strand = ignore.strand,
                    all.a= TRUE)

  # Compute coverage ----
  cov <- idx[, sum(overlap.width>0), keyby= idx.a]$V1

  # Return ----
  return(cov)
}
