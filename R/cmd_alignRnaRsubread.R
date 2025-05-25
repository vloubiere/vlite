#' Generate Commands for Sequence Alignment Using Rsubread
#'
#' @description
#' Creates shell commands to align sequencing reads to a reference genome using the Rsubread aligner.
#' Outputs a BAM file and alignment statistics.
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths.
#' @param fq2 For paired-end data, a character vector of .fq (or .fq.gz) file paths matching fq1 files. Default= NULL.
#' @param output.prefix Prefix for the output files.
#' @param genome Reference genome name (e.g., "mm10", "dm6"). If not provided, genome.idx must be specified.
#' @param genome.idx Path to the Rsubread genome index. Default= NULL.
#' @param bam.output.folder Directory for the output BAM file. Default= "db/bam.output.folder/".
#' @param alignment.stats.output.folder Directory for alignment statistics. Default= "db/alignment_stats/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of CPU cores to use. Default= 8.
#'
#' @return A `data.table` with:
#' - `file.type`: Output file labels ("bam", "align.stats").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run Rsubread.
#' - `cores`: The number of CPU cores to use.
#' - `job.name`: Default name for the job = "alnRsub".
#'
#' @examples
#' # Align single-end reads to the mm10 genome
#' cmd <- cmd_alignRnaRsubread(
#'   fq1 = "sample_R1.fq.gz",
#'   output.prefix = "sample",
#'   genome = "mm10"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Align paired-end reads using a custom genome index
#' cmd <- cmd_alignRnaRsubread(
#'   fq1 = "sample_R1.fq.gz",
#'   fq2 = "sample_R2.fq.gz",
#'   output.prefix = "sample",
#'   genome.idx = "/path/to/custom/genome/index"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_alignRnaRsubread <- function(fq1,
                                 fq2= NULL,
                                 output.prefix,
                                 genome,
                                 genome.idx= NULL,
                                 bam.output.folder= "db/bam.output.folder/",
                                 alignment.stats.output.folder= "db/alignment_stats/",
                                 Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                 cores= 8)
{
  # Check (!Do not check if fq1 or fq2 files exist to allow wrapping!) ----
  fq1 <- unique(fq1)
  if(!is.null(fq2))
    fq2 <- unique(fq2)
  if(any(!grepl(".fq$|.fastq$|.fq.gz$|.fastq.gz$", c(fq1, fq2))))
    stop("fq1 and fq2 file paths should end up with `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`")
  if(!is.null(fq2) && length(fq1) != length(fq2))
    stop("When provided, fq2 files should match fq1 files.")
  if(missing(genome) && is.null(genome.idx))
    stop("genome is missing and and genome.idx is set to NULL -> exit")

  # Retrieve index ----
  if(!missing(genome)) {
    genome.idx <- switch(genome,
                         "mm10"= "/groups/stark/vloubiere/genomes/Mus_musculus/subreadr_mm10/subreadr_mm10_index",
                         "dm6"= "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/subreadr_dm6/subreadr_dm6_index")
  } else {
    genome <- basename(genome.idx)
  }

  # Output files paths ----
  bam <- file.path(bam.output.folder, paste0(output.prefix, "_", genome, ".bam"))
  stats <- paste0(bam, ".summary")

  # Align command ----
  # * If several fq1/fq2 files provided, they will be merged at this step
  cmd <- paste(
    Rpath,
    system.file("Rscript", "align_rna_Rsubread.R", package = "vlite"),
    paste0(fq1, collapse= ","),
    ifelse(is.null(fq2), "''", paste0(fq2, collapse= ",")),
    genome.idx,
    bam
  )

  # Move alignment statistics ----
  stats.new <- file.path(alignment.stats.output.folder, basename(stats))
  cmd <- paste(cmd, "; mv", stats, stats.new)

  # Wrap commands output ----
  cmd <- data.table(file.type= c("bam", "align.stats"),
                    path= c(bam, stats.new),
                    cmd= cmd,
                    cores= cores,
                    job.name= "alnRsub")

  # Return ----
  return(cmd)
}
