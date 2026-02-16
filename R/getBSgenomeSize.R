#' Get BS chromosome sizes
#'
#' @param genome BSgenome object to use: "dm6", "mm10"...
#' @examples
#' getBSchromSizes("mm10")
#'
#' @return data.table containing chromosomes start and end coordinates
#' @export
getBSgenomeSize <- function(genome)
{
  genome <- BSgenome::getBSgenome(genome, load.only = TRUE)
  gSize <- GenomeInfoDb::seqinfo(genome)
  gSize <- data.table(seqnames= GenomicRanges::seqnames(gSize),
                      start= 1L,
                      end= GenomeInfoDb::seqlengths(gSize))
  gSize[, strand:= "*"]
  return(gSize)
}
