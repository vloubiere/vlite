#' bw coverage
#'
#' Compute mean signal from a bigwig file at a set of coordinates. Regions with no overlap will return NA.
#'
#' @param bed BED file input compatible with '?importBed', representing the genomic ranges for which mean bw signal will be returned.
#' @param bw Path to a target bw file.
#'
#' @examples
#' # Example track
#' bw <- "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PH18_merge.bw"
#'
#' Compute coverage
#' bed <- importBed(c("chr2R:10000000-10010000", "chr3R:10000000-10005000", "nonExistingChr:1000-2000"))
#' bwCoverage(bed= bed, bw= bw)
#'
#' @export
bwCoverage <- function(bed,
                       bw)
{
  # Checks ----
  if(length(bw) != 1)
    stop("length(bw) != 1")
  if(!file.exists(bw))
    stop("bw file does not exist!")

  # Copy for incapsulation ----
  bed <- vlite::importBed(bed)

  # Import bw ----
  gr <- GenomicRanges::GRanges(bed)
  sel <- rtracklayer::BigWigSelection(gr, "score")
  var <- suppressWarnings(
    rtracklayer::import.bw(con = bw,
                           selection= sel,
                           as="RleList")
  )
  seqlevels(gr, pruning.mode="coarse") <- names(var)

  # Compute mean score ----
  bed[seqnames %in% seqlevels(gr), cov:= binnedAverage(gr, var, "average_score")$average_score]

  # Return
  return(bed$cov)
}
