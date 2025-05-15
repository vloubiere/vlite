#' Retrieve the seqnames within a bigwig file
#'
#' Returns all the seqLengths of the chromosomes that are represented within a bigwig file.
#'
#' @param bw Path to a target bw file.
#'
#' @examples
#' # Example track
#' bw <- "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PH18_merge.bw"
#' bwGetSeqlengths(bw= bw)
#'
#' @export
bwGetSeqlengths <- function(bw)
{
  # Checks ----
  if(length(bw) != 1)
    stop("length(bw) != 1")
  if(!file.exists(bw))
    stop("bw file does not exist!")

  # Get seqLengths ----
  current <- as.data.frame(GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(track)))
  current <- as.data.table(current, keep.rownames= "seqnames")
  current <- current[, .(seqnames, start= 1, end= seqlengths)]

  # Return
  return(current)
}
