#' Get genomic sequence
#'
#' Returns the sequences of a bed.
#'
#' @param bed Input regions in any format compatible with ?importBed.
#' @param genome BSgenome name: "dm6", "mm10"...
#'
#' @examples
#' # Sample random regions
#' test <- randomRegionsBSgenome("dm6", widths= rep(100, 3))
#' getBSsequence(test, "dm6")
#'
#' @return Character vector of sequences
#' @export
getBSsequence <- function(bed,
                          genome)
{
  # Import ----
  bed <- importBed(bed)

  # Make sequence names ----
  if(!"name" %in% names(bed)) {
    bed[, name:= paste0(seqnames, ":", start, "-", end, ":", strand)]
  }

  # Extract sequences ----
  sequences <- BSgenome::getSeq(BSgenome::getBSgenome(genome, load.only = TRUE),
                                names= bed$seqnames,
                                start= bed$start,
                                end= bed$end,
                                strand= bed$strand,
                                as.character= T)

  # Add names and return ----
  names(sequences) <- bed$name
  return(sequences)
}
