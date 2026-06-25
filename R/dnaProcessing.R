#' Align gDNA using rsubread
#'
#' @description Generates a pipeline of commands to trim gDNA sequencing reads
#' and align them using the rsubread package.
#' 
#' @param fq1 A character vector of .fq (or .fq.gz) file paths.
#' @param fq2 For paired-end data, a character vector of .fq (or .fq.gz) file paths matching fq1 files. Default= NULL.
#' @param output.prefix Prefix for output files.
#' @param genome Reference genome identifier (e.g., "mm10", "hg38").
#' @param genome.idx Path to the Rsubread genome index. If NULL, derived from genome. Default= NULL.
#' @param fq.output.folder Directory for trimmed FASTQ files. Default= "db/fq/DNASeq/".
#' @param bam.output.folder Directory for aligned BAM files. Default= "db/bam/DNASeq/".
#' @param alignment.stats.output.folder Directory for alignment statistics. Default= "db/alignment_stats/DNASeq/".
#' @param Rpath Path to the Rscript binary. Default= "Rscript".
#' @param cores Number of CPU cores to use. Default= 8.
#'
#' @return A data.table with:
#' - `file.types`: Types of output files.
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands for each step in the pipeline.
#' - `cores`: The number of CPU cores to use.
#' - `job.name`: Default name for the job = paste0("DNA_", output.prefix).
#'
#' @examples
#' # Process paired-end RNA-Seq data
#' cmd <- dnaProcessing(
#'   fq1 = c("sample1_R1.fq.gz", "sample1_R1_rep2.fq.gz"),
#'   fq2 = c("sample1_R2.fq.gz", "sample1_R2_rep2.fq.gz"),
#'   output.prefix = "sample1",
#'   genome = "hg38",
#'   cores = 8
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
dnaProcessing <- function(
    fq1,
    fq2= NULL,
    output.prefix,
    genome,
    genome.idx= NULL,
    allowMultiOverlap= FALSE,
    fq.output.folder= "db/fq/DNASeq/",
    bam.output.folder= "db/bam/DNASeq/",
    alignment.stats.output.folder= "db/alignment_stats/DNASeq/",
    counts.stats.output.folder= "db/stats/DNASeq/",
    counts.output.folder= "db/counts/DNASeq/",
    Rpath= "Rscript",
    cores= 8
)
{
  # Trimming illumina adaptors ----
  cmd <- cmd_trimIlluminaAdaptors(
    fq1= fq1,
    fq2= fq2,
    fq.output.folder= fq.output.folder
  )
  
  # * If several fq1/fq2 files provided, they will be merged during alignment ----
  fq1.trim <- cmd[file.type=="fq1.trim", path]
  fq2.trim <- if(!is.null(fq2))
    cmd[file.type=="fq2.trim", path] else
      NULL
  
  # Alignment ----
  align.cmd <- cmd_alignDnaRsubread(
    fq1= fq1.trim,
    fq2= fq2.trim,
    output.prefix= output.prefix,
    genome= genome,
    genome.idx= genome.idx,
    bam.output.folder= bam.output.folder,
    alignment.stats.output.folder = alignment.stats.output.folder,
    Rpath= Rpath,
    cores= cores
  )
  cmd <- rbind(cmd, align.cmd, fill= TRUE)
  
  # # bw tracks ----
  # bw.cmd <- cmd_bamToBigwig(
  #   bam = align.cmd[file.type=="bam", path],
  #   layout = ifelse(is.null(fq2), "SINGLE", "PAIRED"),
  #   output.prefix = NULL, # bam file basename
  #   extend.PE.fragments = FALSE,
  #   libsize.normalize = TRUE,
  #   extsize = 0,
  #   bw.output.folder = bw.output.folder,
  #   Rpath = Rpath,
  #   cores= cores
  # )
  # cmd <- rbind(cmd, bw.cmd, fill= TRUE)
  
  # Return ----
  cmd$cores <- cores
  cmd$job.name <- paste0("DNA_", output.prefix)
  return(cmd)
}
