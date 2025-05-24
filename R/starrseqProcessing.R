#' Process STARR-Seq sequencing Data
#'
#' @description
#' Implements a pipeline for processing ORFeome sequencing data, including:
#' 1. Adapter trimming
#' 2. Genome alignment
#' 3. UMI collapsing
#' 4. Bigwig tracks
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths.
#' @param fq2 For paired-end data, a character vector of .fq (or .fq.gz) file paths matching fq1 files. Default= NULL.
#' @param output.prefix Prefix for output files.
#' @param genome Reference genome identifier (e.g., "mm10", "hg38").
#' @param genome.idx Path to the Rsubread genome index. If NULL, derived from genome. Default= NULL.
#' @param max.ins Bowtie2 maximum insert size for paired-end alignment. Default= 2000.
#' @param bed.subset Optional path to a bed file. If provided, only reads overlapping these regions on the same strand are used to generate
#' bigwig tracks. Default= NULL (all reads are used).
#' @param fq.output.folder Directory where trimmed fastq files will be saved. Default= "db/fq/STARRSeq/".
#' @param bam.output.folder Directory where aligned bam files will be saved. Default= "db/bam/STARRSeq/".
#' @param umi.counts.output.folder Directory where UMI counts will be saved. Default= "db/umi_counts/STARRSeq/".
#' @param bed.output.folder Directory where UMI-collapsed bed files will be saved. Default= "db/bed/STARRSeq/".
#' @param bw.output.folder Directory where bigwig tracks will be saved. Default= "db/bw/STARRSeq/".
#' @param Rpath Path to the Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of CPU cores to use. Default= 8.
#'
#' @return A data.table with:
#' - `file.type`: Type of output file.
#' - `path`: Path to the output file.
#' - `cmd`: Shell command for each step in the pipeline.
#' - `cores`: The number of CPU cores to use.
#' - `job.name`: Default name for the job = paste0("STARR_", output.prefix).
#'
#' @examples
#' # Process ORFeome sequencing data
#' cmd <- starrseqProcessing(
#'   fq1 = c("sample1_R1.fq.gz", "sample1_R1_rep2.fq.gz"),
#'   output.prefix = "sample1",
#'   cores = 8
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Then, peaks should be called using ?cmd_MAGeCK_ORFeome()
#'
#' @export
starrseqProcessing <- function(fq1,
                               fq2= NULL,
                               output.prefix,
                               genome,
                               genome.idx= NULL,
                               max.ins= 2000,
                               bed.subset= NULL,
                               fq.output.folder= "db/fq/STARRSeq/",
                               bam.output.folder= "db/bam/STARRSeq/",
                               alignment.stats.output.folder= "db/alignment_stats/STARRSeq/",
                               umi.counts.output.folder= "db/umi_counts/STARRSeq/",
                               bed.output.folder= "db/bed/STARRSeq/",
                               bw.output.folder= "db/bw/STARRSeq/",
                               Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                               cores= 8)
{
  # Trimming illumina adaptors ----
  cmd <- cmd_trimIlluminaAdaptors(fq1 = fq1,
                                  fq2 = fq2,
                                  fq.output.folder = fq.output.folder)

  # * If several fq1/fq2 files provided, they will be merged during alignment ----
  fq1.trim <- cmd[file.type=="fq1.trim", path]
  fq2.trim <- if(!is.null(fq2))
    cmd[file.type=="fq2.trim", path] else
      NULL

  # Alignment ----
  align.cmd <- cmd_alignBowtie2(fq1 = fq1.trim,
                                fq2= fq2.trim,
                                output.prefix = output.prefix,
                                genome = genome,
                                genome.idx = genome.idx,
                                max.ins = max.ins,
                                mapq = 30,
                                bam.output.folder = bam.output.folder,
                                alignment.stats.output.folder = alignment.stats.output.folder,
                                cores = cores)
  cmd <- rbind(cmd, align.cmd, fill= TRUE)

  # UMI counts ----
  umi.cmd <- cmd_umiCountsFromBam(bam = cmd[file.type=="bam", path],
                                  layout= ifelse(is.null(fq2), "SINGLE", "PAIRED"),
                                  output.prefix= output.prefix,
                                  umi.counts.output.folder= umi.counts.output.folder,
                                  collapsed.bed.output.folder= bed.output.folder)
  cmd <- rbind(cmd, umi.cmd, fill= TRUE)

  # Generate bw tracks ----
  bw.cmd <- cmd_bedToBigwig(bed = cmd[file.type=="umi.bed", path],
                            genome = genome,
                            bed.subset = bed.subset,
                            output.prefix = output.prefix,
                            bw.output.folder = bw.output.folder)
  cmd <- rbind(cmd, bw.cmd, fill= TRUE)

  # Return ----
  cmd$cores <- cores
  cmd$job.name <- paste0("STARR_", output.prefix)
  return(cmd)
}
