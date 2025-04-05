#' Get BS genome sizes
#'
#' Very fast using data.table
#'
#' @param genome BSgenome object to use. Example "dm3", "mm10"...
#' @examples
#' getBSgenomeSize("mm10")
#'
#' @return data.table bed containing chromosomes start and (seqLength) coordinates
#' @export
getBSgenomeSize <- function(genome)
{
  gSize <- GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome, load.only = TRUE))
  gSize <- data.table(seqnames= seqnames(gSize),
                      start= 1L,
                      end= seqlengths(gSize))
  return(gSize)
}
