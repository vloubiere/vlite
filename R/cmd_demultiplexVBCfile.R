#' Generate Demultiplexing Commands for .bam or .tar.gz Files from the VBC NGS Facility
#'
#' @description
#' Generates shell commands to demultiplex BAM or tar.gz FASTQ files from the VBC NGS facility.
#' This function is a wrapper around the Perl scripts vbc_bam_demultiplexing.pl and vbc_tar_demultiplexing.pl.
#' It supports both single-end and paired-end sequencing data, with optional support for PRO-Seq-specific processing
#' (only supported with BAM input).
#'
#' @param vbcFile Path to the input BAM file or .tar.gz file containing the reads.
#' @param layout Sequencing layout, must be either "PAIRED" or "SINGLE".
#' @param i7 i7 index sequence(s). Multiple indexes should be comma-separated.
#'        If set to "none" (default), no i7 filtering is performed
#' @param i5 i5 index sequence(s). Multiple indexes should be comma-separated.
#'        If set to "none" (default), no i5 filtering is performed
#' @param i7.column Column number in the BAM file containing the i7 index. Default= 14L.
#' @param i5.column Column number in the BAM file containing the i5 index. Default= 12L.
#' @param output.prefix Output files prefix. If not provided,
#'        constructed from the input file name and index sequences.
#' @param start.seq For PRO-Seq data, the eBC DNA string that must be present at the start of the reads.
#'        Only supported for BAM input. Default= NULL.
#' @param trim.length When start.seq is provided (PRO-Seq reads), number of nucleotides that should be cut
#'        from the sequence (after start.seq has been trimmed) and appended to the read ID. Only supported for BAM input.
#' @param fq.output.folder Directory where output FASTQ files will be written. Default= "db/fq/".
#' @param cores Number of CPU cores to use for samtools processing (when using BAM input). Default= 8L.
#' @param head Number of reads that should be processed (for testing purposes).
#'
#' @return A data.table with three columns:
#' \itemize{
#'   \item `file.type`: Labels for output files (e.g.: "fq1", "fq2").
#'   \item `path`: Full paths to the output files.
#'   \item `cmd`: Shell command(s) to generate the output files.
#' }
#'
#' @details
#' The function generates commands that:
#' 1. Read a BAM file using samtools or a .tar.gz file.
#' 2. Process the reads using vbc_bam_demultiplexing.pl (BAM) or vbc_tar_demultiplexing.pl (tar.gz).
#' 3. Filter reads based on i7 and/or i5 indexes.
#' 4. Optionally process PRO-Seq-specific requirements (BAM only).
#' 5. Output gzipped FASTQ files.
#'
#' @section Output Files:
#' For paired-end data:
#' - `<output.prefix>_1.fq.gz`
#' - `<output.prefix>_2.fq.gz`
#'
#' For single-end data:
#' - `<output.prefix>.fq.gz`
#'
#' @examples
#' # Example with a .tar.gz file (One of Tomas screens)
#' cmd <- cmd_demultiplexVBCfile(
#'   vbcFile = "/groups/stark/projects/PE75_20250120_TP/2415402230_1_R18342_20250120.tar.gz",
#'   layout = "PAIRED",
#'   i7 = "CTATAC",
#'   i5 = "none",
#'   output.prefix = "ORFeome_fastq_demult_CTATAC",
#'   fq.output.folder = "/groups/stark/vloubiere/packages/tests/",
#'   head = 40000
#' )
#' vl_submit(cmd, overwrite = FALSE, logs = "/groups/stark/vloubiere/packages/tests/logs/")
#'
#' # Example with a BAM file (PRO-Seq data)
#' cmd <- cmd_demultiplexVBCfile(
#'   vbcFile = "/groups/stark/projects/PE50_20230401/AAAYLY5HV_1_20230331B_20230401.bam",
#'   layout = "PAIRED,
#'   i7 = "GTCCGC",
#'   i5 = "none,
#'   i7.column = 14,
#'   output.prefix = "PE_PROSeq",
#'   fq.output.folder = "/groups/stark/vloubiere/packages/tests/",
#'   start.seq = "ATCG",
#'   trim.length = 10,
#'   cores = 8,
#'   head = 40000
#' )
#' vl_submit(cmd,
#'           overwrite = TRUE,
#'           logs = "/groups/stark/vloubiere/packages/tests/logs/")
#'
#' @seealso
#' The underlying Perl script documentation can be accessed using:
#' \code{system2(system.file("demultiplexing", "vbc_bam_demultiplexing.pl",
#'              package = "genomicsPipelines"), args = "--help")}
#'
#' @export
cmd_demultiplexVBCfile <- function(vbcFile,
                                   layout,
                                   i7= "none",
                                   i5= "none",
                                   umi= FALSE,
                                   i7.column= 14,
                                   i5.column= 12,
                                   output.prefix,
                                   fq.output.folder= "db/fq/",
                                   start.seq= NULL,
                                   trim.length,
                                   cores= 8,
                                   head)
{
  options(scipen= 999)
  # Checks ----
  if(length(vbcFile)!=1)
    stop("A unique vbcFile path should be provided.")
  if(!file.exists(vbcFile))
    stop("vbcFile could not be found.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(length(i7)>1 | length(i5)>1)
    stop("If several i7 or i5 indexes have been used for one sample, they should be concatenated and comma-separated.")
  if(length(i7.column)>1 | length(i5.column)>1)
    stop("Several i7 or i5 column number were provided.")
  if(!is.null(start.seq) && !grepl(".bam$", vbcFile))
    stop("start.seq is only supported for bam files.")
  if(!is.null(start.seq) && missing(trim.length))
    stop("When start.seq is specified (PRO-Seq reads), trim.length should also be provided.")
  if(!missing(head) && head %% 1!=0)
    stop("head should be a round number.")

  # Default output prefix ----
  if(missing(output.prefix)) {
    output.prefix <- paste0(gsub(".bam$|.tar.gz$", "", basename(vbcFile)),
                            "_", gsub(",", ".", i7), "_", gsub(",", ".", i5))
  }

  # Output files paths ----
  output.prefix <- file.path(fq.output.folder, output.prefix)
  fq1 <- paste0(output.prefix, ifelse(layout=="PAIRED", "_1.fq.gz", ".fq.gz"))
  if(layout=="PAIRED")
    fq2 <- paste0(output.prefix, "_2.fq.gz")

  # Method for bam files ----
  if(grepl(".bam$", vbcFile)) {
    # Decompress command
    decompress <- paste("samtools view -@", cores-1, vbcFile, "|")
    if(!missing(head))
      decompress <- paste(decompress, "head -n", head, "|")

    # Demultiplexing command
    cmd <- paste(
      decompress,
      "perl", system.file("perl", "vbc_bam_demultiplexing.pl", package = "genomicsPipelines"), # perl script
      paste0("'", layout, "'"),
      paste0("'", i7, "'"), # i7 barcode sequence
      paste0("'", i5, "'"), # i5 index sequence
      i7.column, # i7 column (default= 14)
      i5.column, # i5 column (default= 12)
      paste0("'", output.prefix, "'")
    )

    # PRO-Seq specific arguments
    if(!is.null(start.seq)) {
      cmd <- paste(cmd,
                   paste0("'", start.seq, "'"),
                   trim.length)
    }
  } else if(grepl(".tar.gz$", vbcFile)) {

    # Method for .tar.gz files ----

    # perl script call
    perl.script <- system.file("perl", "vbc_tar_demultiplexing.pl", package = "genomicsPipelines")
    perl.script <- paste("perl", perl.script)
    if(!missing(head))
      perl.script <- paste(perl.script, "-n", head)

    # Demultiplexing parameters
    cmd <- paste(perl.script,
                 ifelse(umi, "--umi", ""),
                 vbcFile,
                 paste0("'", i7, "'"), # i7 barcode sequence
                 paste0("'", i5, "'"), # i5 index sequence
                 paste0("'", output.prefix, "'"),
                 paste0("'", layout, "'"))

    # Check PRO-Seq args
    if(!is.null(start.seq))
      stop("Optional'start.seq' parameter for PRO-Seq is not supoorted for .tar.gz input files.")
  } else {
    stop("vbcFile should be in .bam or .tar.gz format.")
  }

  # Wrap commands output ----
  cmd <- data.table(file.type= if(layout=="PAIRED") c("fq1", "fq2") else "fq1",
                    path= if(layout=="PAIRED") c(fq1, fq2) else fq1,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
