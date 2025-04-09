#' Process ORFtag Sequencing Data
#'
#' @description
#' Implements a pipeline for processing ORFtag sequencing data, including:
#' 1. Adapter trimming
#' 2. Genome alignment
#' 3. BAM file collapsing (unique insertions)
#' 4. Assigning insertions to genomic features
#'
#' Supports single-end sequencing data, with options for custom genome and annotation files.
#'
#' @param fq1 Path(s) to input R1 FASTQ file(s). Note that R2 is not used here.
#' @param output.prefix Prefix for output files.
#' @param genome Reference genome identifier (e.g., `"mm10"`, `"hg38"`).
#' @param genome.idx Path to the Bowtie2 genome index. If `NULL`, derived from `genome`. Default: `NULL`.
#' @param gtf Path to the GTF annotation file. Default: `NULL`.
#' @param fq.output.folder Directory for trimmed FASTQ files. Default: `"db/fq/ORFtag/"`.
#' @param bam.output.folder Directory for aligned BAM files. Default: `"db/bam/ORFtag/"`.
#' @param alignment.stats.output.folder Directory for alignment statistics. Default: `"db/alignment_stats/ORFtag/"`.
#' @param bed.output.folder Directory for BED files of unique insertions. Default: `"db/bed/ORFtag/"`.
#' @param counts.output.folder Directory for counts files. Default: `"db/counts/ORFtag/"`.
#' @param Rpath Path to the Rscript binary. Default: `"/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"`.
#' @param cores Number of CPU cores to use. Default: `8`.
#'
#' @return A `data.table` with:
#' - `file.types`: Types of output files.
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands for each step in the pipeline.
#'
#' @examples
#' # Process ORFtag sequencing data
#' cmd <- orftagProcessing(
#'   fq1 = c("sample1_R1.fq.gz", "sample1_R1_rep2.fq.gz"),
#'   output.prefix = "sample1",
#'   genome = "mm10",
#'   gtf = "/data/annotations/mm10.gtf",
#'   cores = 8
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{cmd_trimIlluminaAdaptors}} for adapter trimming
#'   \item \code{\link{cmd_alignBowtie2}} for alignment
#'   \item \code{\link{cmd_collapseBam}} for BAM file collapsing
#'   \item \code{\link{cmd_assignInsertions}} for assigning insertions
#' }
#'
#' @export
orftagProcessing <- function(fq1,
                             output.prefix,
                             genome,
                             genome.idx= NULL,
                             gtf= NULL,
                             fq.output.folder= "db/fq/ORFtag/",
                             bam.output.folder= "db/bam/ORFtag/",
                             alignment.stats.output.folder= "db/alignment_stats/ORFtag/",
                             bed.output.folder= "db/bed/ORFtag/",
                             counts.output.folder= "db/counts/ORFtag/",
                             Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                             cores= 8)
{
  # Initiate command output ----
  cmd <- data.table(file.type= character(),
                    path= character(),
                    cmd= character())

  # Trimming illumina adaptors ----
  for(i in seq(fq1)) { # Multiple files can be provided for one sample
    .c <- cmd_trimIlluminaAdaptors(fq1= fq1[i],
                                   fq2= NULL, # Second reads are not used
                                   fq.output.folder= fq.output.folder)
    cmd <- rbind(cmd, .c)
  }

  # * If several fq1 files provided, they will be merged at this step ----
  fq1.trim <- paste0(cmd[file.type=="fq1.trim", path], collapse = ",")

  # Alignment ----
  align.cmd <- cmd_alignBowtie2(fq1= fq1.trim,
                                fq2= NULL, # Second reads are not used
                                output.prefix= output.prefix,
                                genome= genome,
                                genome.idx= genome.idx,
                                mapq= 30,
                                bam.output.folder= bam.output.folder,
                                alignment.stats.output.folder= alignment.stats.output.folder,
                                cores= cores)
  cmd <- rbind(cmd, align.cmd)

  # Collapse bam file (unique insertions) ----
  collapse.cmd <- cmd_collapseBam(bam = cmd[file.type=="bam", path],
                                  output.prefix = NULL, # From bam file
                                  collapsed.bam.output.folder = bam.output.folder,
                                  collapsed.stats.output.folder = alignment.stats.output.folder,
                                  cores = cores)
  cmd <- rbind(cmd, collapse.cmd)

  # Assign insertions to closest downstream exon ----
  cmd <- cmd_assignInsertions(bam = cmd[file.type=="collapsed.bam", path],
                              output.prefix = NULL, # From bam file
                              genome = genome,
                              gtf = gtf,
                              bed.output.folder = bed.output.folder,
                              counts.output.folder = counts.output.folder,
                              Rpath = Rpath)
  cmd <- rbind(cmd, collapse.cmd)

  # Return ----
  return(cmd)
}
