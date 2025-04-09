#' Generate Cutadapt Commands for Trimming PRO-Seq Adapters
#'
#' @description
#' Creates shell commands to trim PRO-Seq adapters from gzipped FASTQ files using Cutadapt.
#'
#' @param fq1 Path to the input FASTQ file (.fq.gz). Note that only read 1 is used in PROseq.
#' @param fq.output.folder Directory for trimmed FASTQ files. Default: `"db/fq/"`.
#'
#' @return A `data.table` with:
#' - `file.type`: Output file label.
#' - `path`: Path to the trimmed FASTQ file.
#' - `cmd`: Shell command to run Cutadapt.
#'
#' @examples
#' cmd <- cmd_trimProseqAdaptors(fq1 = "sample_R1.fq.gz", fq.output.folder = "db/fq/")
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_trimProseqAdaptors <- function(fq1,
                                   fq.output.folder= "db/fq/")
{
  # Check ----
  if(length(fq1)!=1)
    stop("A unique fq1 file should be provided.")

  # Output trimmed fq files paths ----
  fq1.trim <- file.path(fq.output.folder, gsub(".fq.gz", "_trimmed.fq.gz", basename(fq1)))

  # Create command ----
  cmd <- paste("cutadapt -a TGGAATTCTCGGGTGCCAAGG",
               "-m", 10, # Remove trimmed reads shorter than 10bp
               "-o", fq1.trim, fq1)

  # Wrap commands output ----
  cmd <- data.table(file.type= "fq1.trim",
                    path= fq1.trim,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
