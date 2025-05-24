#' Generate Cutadapt Commands for Trimming PRO-Seq Adapters
#'
#' @description
#' Creates shell commands to trim PRO-Seq adapters from gzipped FASTQ files using Cutadapt.
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths.
#' @param fq.output.folder Directory for trimmed FASTQ files. Default= "db/fq/".
#' @param min.length Remove reads shorter than trim.lengt (after trimming). Defaul= 10.
#'
#' @return A `data.table` with:
#' - `file.type`: Output file label.
#' - `path`: Path to the trimmed FASTQ file.
#' - `cmd`: Shell command to run Cutadapt.
#' - `job.name`: Default name for the job = "trimCut".
#'
#' @examples
#' cmd <- cmd_trimProseqAdaptors(fq1 = "sample_R1.fq.gz", fq.output.folder = "db/fq/")
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_trimProseqAdaptors <- function(fq1,
                                   fq.output.folder= "db/fq/",
                                   min.length= 10)
{
  # Check (!Do not check if fq1 or fq2 files exist to allow wrapping!) ----
  fq1 <- unique(fq1)
  if(!grepl(".fq$|.fastq$|.fq.gz$|.fastq.gz$", fq1))
    stop("fq1 file paths should end up with `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`")

  # Output trimmed fq files paths ----
  fq1.trim <- file.path(fq.output.folder,
                        paste0(sub("\\.(fq|fastq)(\\.gz)?$", "", basename(fq1)), "_trimmed.fq.gz"))

  # Create command ----
  cmd <- sapply(seq(fq1), function(i) {
    paste("cutadapt -a TGGAATTCTCGGGTGCCAAGG",
          "-m", min.length, # Remove trimmed reads shorter than 10bp
          "-o", fq1.trim[i], # Output file
          fq1[i]) # Input file
  })

  # Wrap commands output ----
  cmd <- data.table(file.type= "fq1.trim",
                    path= fq1.trim,
                    cmd= cmd,
                    job.name= "trimCut")

  # Return ----
  return(cmd)
}
