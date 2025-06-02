#' motifPosToMatrix
#'
#' Transforms the output of vl_motifPos() to a list of scores' matrices.
#'
#' @param motifPos A ?vl_motifPos() output data.table.
#' @param seqWidth. The width of the sequences for which motif positions have been computed.
#' @param cleanup.cache If set to TRUE, overwrites cached intermediate results. Default= FALSE.
#'
#' @return Returns a list of maximum scores' matrices (one per motif).
#' @export
#'
#' @examples
motifPosToBed <- function(motifPos)
{
  # Checks
  if(!is.data.table(motifPos))
    stop("motifPos should be a data.table similar to ?vl_motifPos() output.")
  if(!all(c("motif", "seqlvls", "ir") %in% names(motifPos)))
    stop("motifPos should contain 'motif', 'seqlvls', 'mot.count' and 'ir' columns (see ?vl_motifPos()).")

  # Melt ----
  bed <- motifPos[mot.count>0, ir[[1]], .(seqlvls, motif)]

  # Return
  return(bed)
}
