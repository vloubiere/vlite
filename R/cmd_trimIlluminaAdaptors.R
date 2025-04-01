#' Generate Commands for Trimming Illumina Adapters Using Trim Galore
#'
#' @description
#' Generates shell commands to trim Illumina adapters from FASTQ files using Trim Galore.
#' This function processes single FASTQ files or pairs of FASTQ files for both single-end
#' and paired-end sequencing data.
#'
#' @param fq1 character(1). Path to the input FASTQ file for single-end data,
#'        or the first read file for paired-end data. Must be gzipped (.fq.gz).
#' @param fq2 character(1). Path to the second read file for paired-end data.
#'        Must be gzipped (.fq.gz). Default: NULL.
#' @param fq.output.folder character(1). Directory where trimmed FASTQ files will be written.
#'
#' @return A data.table with three columns:
#' \itemize{
#'   \item file.type: Labels for output files ("fq1.trim", "fq2.trim" for paired-end,
#'         or "fq1" for single-end)
#'   \item path: Full paths to the output files
#'   \item cmd: Shell command to run Trim Galore
#' }
#'
#' @details
#' The function generates commands that will:
#' 1. Process input FASTQ file(s) using Trim Galore
#' 2. Remove Illumina adapter sequences
#' 3. Output gzipped, trimmed FASTQ files
#'
#' For paired-end data, the output files will be named:
#' - <fq1_basename>_val_1.fq.gz
#' - <fq2_basename>_val_2.fq.gz
#'
#' For single-end data, the output file will be named:
#' - <fq1_basename>_trimmed.fq.gz
#'
#' @section Requirements:
#' - Trim Galore must be installed and available in the system PATH
#' - Input FASTQ files must be gzipped (.fq.gz)
#' - For paired-end data, both files must exist and be properly paired
#'
#' @examples
#' \dontrun{
#' # Single-end processing
#' cmd_trimIlluminaAdaptors(
#'   fq1 = "sample_R1.fq.gz",
#'   fq.output.folder = "trimmed_output"
#' )
#'
#' # Paired-end processing
#' cmd_trimIlluminaAdaptors(
#'   fq1 = "sample_R1.fq.gz",
#'   fq2 = "sample_R2.fq.gz",
#'   fq.output.folder = "trimmed_output"
#' )
#' }
#'
#' @seealso
#' Trim Galore documentation:
#' \url{https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/}
#'
#' @section Output Files:
#' The function generates trimmed FASTQ files with the following naming convention:
#' \itemize{
#'   \item Paired-end: *_val_1.fq.gz and *_val_2.fq.gz
#'   \item Single-end: *_trimmed.fq.gz
#' }
#'
#' @section Warning:
#' This function is not parallelized and is designed to process a single input file
#' or a single pair of FASTQ files at a time.
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
  if(!dir.exists(fq.output.folder))
    dir.create(fq.output.folder, recursive = TRUE, showWarnings = FALSE)

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
