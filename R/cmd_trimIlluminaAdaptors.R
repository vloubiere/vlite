#' Generate Trim Galore Commands for FASTQ Files
#'
#' @description
#' Creates shell commands to trim Illumina adapters from gzipped FASTQ files using Trim Galore.
#' Supports single-end and paired-end sequencing data.
#'
#' @param fq1 Path to the input FASTQ file (.fq.gz) for single-end data or the first read file for paired-end data.
#' @param fq2 Path to the second read file for paired-end data (.fq.gz). Default: `NULL`.
#' @param fq.output.folder Directory for trimmed FASTQ files. Default: `"db/fq/"`.
#'
#' @return A `data.table` with:
#' - `file.type`: Output file labels.
#' - `path`: Paths to trimmed FASTQ files.
#' - `cmd`: Shell command to run Trim Galore.
#'
#' @examples
#' # Single-end
#' cmd <- cmd_trimIlluminaAdaptors(fq1 = "sample_R1.fq.gz")
#' vl_submit(cmd, execute= FALSE)
#'
#' # Paired-end
#' cmd <- cmd_trimIlluminaAdaptors(fq1 = "sample_R1.fq.gz", fq2 = "sample_R2.fq.gz")
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_trimIlluminaAdaptors <- function(fq1,
                                     fq2= NULL,
                                     fq.output.folder= "db/fq/")
{
  # Check ----
  if(length(fq1)!=1)
    stop("A unique fq1 file should be provided.")
  if(!is.null(fq2) && length(fq2)!=1)
    stop("A unique fq2 file should be provided.")

  # Output trimmed fq files paths and commands ----
  if(!is.null(fq2)) {
    # Paired-end reads
    cmd <- paste("trim_galore --gzip --paired -o", fq.output.folder, fq1, fq2)
    fq1.trim <- file.path(fq.output.folder, gsub("(.*)(.f.*q.*)", "\\1_val_1\\2", basename(fq1)))
    fq2.trim <- file.path(fq.output.folder, gsub("(.*)(.f.*q.*)", "\\1_val_2\\2", basename(fq2)))
  } else {
    # Single-end reads
    cmd <- paste("trim_galore --gzip -o", fq.output.folder, fq1)
    fq1.trim <- file.path(fq.output.folder, gsub("(.*)(.f.*q.*)", "\\1_trimmed\\2", basename(fq1)))
  }

  # Wrap commands output ----
  cmd <- data.table(file.type= if(!is.null(fq2)) c("fq1.trim", "fq2.trim") else "fq1.trim",
                    path= if(!is.null(fq2)) c(fq1.trim, fq2.trim) else fq1.trim,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
