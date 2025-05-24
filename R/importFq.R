#' importFq
#'
#' A wrapper to import fq files using data.table::fread().
#'
#' @param fq.path Path to the .fq of .fq.gz file to import.
#' @param head An optional integer specifying the number of reads that should be imported.
#'
#' @return
#' A data.table with columns:
#'  - readID
#'  - Sequence
#' @export
#'
#' @examples
importFq <- function(fq.path, head)
{
  if(length(fq.path)!=1)
    stop("fq.path should be unique.")
  if(!grepl(".fq$|.fastq$|.fq.gz$|.fastq.gz$", fq.path))
    stop("fq.path should be in format .fq, .fastq, .fq.gz, or .fastq.gz.")

  # Pipe the file
  cmd <- if(grepl(".gz$", fq.path))
    paste("zcat", fq.path, "|") else
      paste(fq.path, "|")
  # Head
  if(!missing(head))
    cmd <- paste(cmd , "head -n", head*4, "|")
  # Parse command
  cmd <- paste(cmd, "perl", system.file("perl", "parseFq.pl", package = "vlite"))
  # Import
  fread(cmd= cmd)
}
