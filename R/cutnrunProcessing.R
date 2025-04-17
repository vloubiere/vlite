#' Process CUT&RUN Sequencing Data
#'
#' @description
#' Implements a pipeline for processing CUT&RUN sequencing data, including:
#' 1. Adapter trimming
#' 2. Genome alignment and quality filtering
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths.
#' @param fq2 For paired-end data, a character vector of .fq (or .fq.gz) file paths matching fq1 files. Default= NULL.
#' @param output.prefix Prefix for output files.
#' @param genome Reference genome identifier (e.g., "mm10", "hg38").
#' @param genome.idx Path to Bowtie2 index. If NULL, derived from genome. Default= NULL.
#' @param fq.output.folder Directory for trimmed FASTQ files. Default= "db/fq/CUTNRUN/".
#' @param bam.output.folder Directory for aligned BAM files. Default= "db/bam/CUTNRUN/".
#' @param alignment.stats.output.folder Directory for alignment statistics. Default= "db/alignment_stats/CUTNRUN/".
#' @param cores Number of CPU cores to use. Default= 8.
#'
#' @return A data.table with:
#' - `file.type`: Types of output files.
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands for each step in the pipeline.
#'
#' @examples
#' # Process paired-end data from FASTQ files
#' cmd <- cutnrunProcessing(
#'   fq1 = c("sample1_R1.fq.gz", "sample1_R1_rep2.fq.gz"),
#'   fq2 = c("sample1_R2.fq.gz", "sample1_R2_rep2.fq.gz"),
#'   output.prefix = "sample1",
#'   genome = "hg38",
#'   cores = 8
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Then, peaks should be called using ?cmd_peakCalling() and ?cmd_confidentPeaks()
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{cmd_trimIlluminaAdaptors}} for adapter trimming
#'   \item \code{\link{cmd_alignBowtie2}} for alignment
#' }
#'
#' @export
cutnrunProcessing <- function(fq1,
                              fq2= NULL,
                              output.prefix,
                              genome,
                              genome.idx= NULL,
                              fq.output.folder= "db/fq/CUTNRUN/",
                              bam.output.folder= "db/bam/CUTNRUN/",
                              alignment.stats.output.folder= "db/alignment_stats/CUTNRUN/",
                              cores= 8)
{
  # Trimming illumina adaptors ----
  cmd <- cmd_trimIlluminaAdaptors(fq1= fq1,
                                  fq2= fq2,
                                  fq.output.folder= fq.output.folder)

  # * If several fq1/fq2 files provided, they will be merged during alignment ----
  fq1.trim <- cmd[file.type=="fq1.trim", path]
  if(!is.null(fq2))
    fq2.trim <- cmd[file.type=="fq2.trim", path]

  # Alignment ----
  align.cmd <- cmd_alignBowtie2(fq1= fq1.trim,
                                fq2= fq2.trim,
                                output.prefix= output.prefix,
                                genome= genome,
                                genome.idx= genome.idx,
                                mapq= 30,
                                bam.output.folder= bam.output.folder,
                                alignment.stats.output.folder= alignment.stats.output.folder,
                                cores= cores)
  cmd <- rbind(cmd, align.cmd)

  # Return ----
  return(cmd)
}
