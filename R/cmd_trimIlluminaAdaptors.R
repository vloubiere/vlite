#' Generate Trim Galore Commands for FASTQ Files
#'
#' @description
#' Creates shell commands to trim Illumina adapters from gzipped FASTQ files using Trim Galore.
#' Supports single-end and paired-end sequencing data.
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths.
#' @param fq2 For paired-end data, a character vector of .fq (or .fq.gz) file paths matching fq1 files. Default= NULL.
#' @param cores The number of cores to use.
#' @param fq.output.folder Directory for trimmed FASTQ files. Default= "db/fq/".
#'
#' @return A `data.table` with:
#' - `file.type`: Output file labels.
#' - `path`: Paths to trimmed FASTQ files.
#' - `cmd`: Shell command to run Trim Galore.
#' - `job.name`: Default name for the job = "trimGal".
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
                                     cores= 1,
                                     fq.output.folder= "db/fq/")
{
  # Check (!Do not check if fq1 or fq2 files exist to allow wrapping!) ----
  fq1 <- unique(fq1)
  if(!is.null(fq2))
    fq2 <- unique(fq2)
  if(any(!grepl(".fq$|.fastq$|.fq.gz$|.fastq.gz$", c(fq1, fq2))))
    stop("fq1 and fq2 file paths should end up with `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`")
  if(!is.null(fq2) && length(fq1) != length(fq2))
    stop("When provided, fq2 files should match fq1 files.")

  # Output file paths ----
  if(is.null(fq2)) {
    # Single-end reads
    fq1.trim <- file.path(fq.output.folder,
                          paste0(sub("\\.(fq|fastq)(\\.gz)?$", "", basename(fq1)), "_trimmed.fq.gz"))
  } else {
    # Paired-end reads
    fq1.trim <- file.path(
      fq.output.folder,
      paste0(sub("\\.(fq|fastq)(\\.gz)?$", "", basename(fq1)), "_val_1.fq.gz")
    )
    fq2.trim <- file.path(
      fq.output.folder,
      paste0(sub("\\.(fq|fastq)(\\.gz)?$", "", basename(fq2)), "_val_2.fq.gz")
    )
  }

  # Trimming commands ----
  cmd <- if(is.null(fq2)) {
    # Single-end reads
    sapply(fq1, function(x) paste("trim_galore --cores", cores, "--gzip -o", fq.output.folder, x))
  } else {
    # Paired-end reads
    sapply(seq(fq1), function(i) paste("trim_galore --cores", cores, "--gzip --paired -o", fq.output.folder, fq1[i], fq2[i]))
  }

  # Wrap commands output ----
  cmd <- data.table(
    file.type= if(is.null(fq2)) "fq1.trim" else rep(c("fq1.trim", "fq2.trim"), each= length(fq1)),
    path= if(is.null(fq2)) fq1.trim else c(fq1.trim, fq2.trim),
    cmd= cmd,
    job.name= "trimGal"
  )

  # Return ----
  return(cmd)
}
