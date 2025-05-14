#' Retrieve the seqnames within a bigwig file
#'
#' Returns all the seqnames that are represented within a bed file.
#'
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
bwGetSeqnames <- function(bw)
{
  # Checks ----
  if(length(bw) != 1)
    stop("length(bw) != 1")
  if(!file.exists(bw))
    stop("bw file does not exist!")

  # Get chromosomes ----
  chr <- seqnames(GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(track)))

  # Return
  return(chr)
}
