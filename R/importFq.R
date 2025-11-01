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
importFq <- function(fq.path, head= NULL)
{
  if(length(fq.path)!=1)
    stop("fq.path should be unique.")
  if(!grepl(".fq$|.fastq$|.fq.gz$|.fastq.gz$", fq.path))
    stop("fq.path should be in format .fq, .fastq, .fq.gz, or .fastq.gz.")

  # Import
  dat <- if(!is.null(head))
    fread(fq.path, header= F, nrows = head) else
      fread(fq.path, header= F)

  # Format
  dat <- data.table(read_ID= dat[seq(.N) %% 4 == 1, V1],
                    seq= dat[seq(.N) %% 4 == 2, V1])

  # Return
  return(dat)
}
