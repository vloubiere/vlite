#' Get genomic sequence
#'
#' Returns the sequences of a bed.
#'
#' @param bed Bed file for which regions have to be returned
#' @param genome BSgenome ID ("dm3", "dm6"...)
#'
#' @examples
#' getBSsequence(vl_SUHW_top_peaks, "dm3")
#'
#' @return Character vector of sequences
#' @export
getBSsequence <- function(bed,
                          genome)
{
  # Import ----
  bed <- importBed(bed)

  # Checks ----
  if(!"end" %in% names(bed)){
    bed[, end:= start]
    warnings("Input 'bed' did not have an 'end' column, which was set to start.")
  }
  if(!"strand" %in% names(bed)){
    bed[, strand:= "*"]
    warnings("Input 'bed' did not have a 'strand' column, which was set to '*'.")
  }

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
