#' Process ORFeome Sequencing Data
#'
#' @description
#' Implements a pipeline for processing ORFeome sequencing data, including:
#' 1. Adapter trimming
#' 2. Genome alignment
#' 3. Barcode counting
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths. Note that read 2 is not used.
#' @param output.prefix Prefix for output files.
#' @param bowtie2.lib.idx The bowtie2 barcodes index. Default= "/groups/stark/vloubiere/projects/viralORF_tomas/db/indexes_BCs/lib200_merged/ORF".
#' @param fq.output.folder Directory for trimmed FASTQ files. Default= "db/fq/ORFeome/".
#' @param bam.output.folder Directory for aligned BAM files. Default= "db/bam/ORFeome/".
#' @param counts.output.folder Directory for barcode counts. Default= "db/counts/ORFeome/".
#' @param alignment.stats.output.folder Directory for alignment statistics. Default= "db/alignment_stats/ORFeome/".
#' @param Rpath Path to the Rscript binary. Default= "Rscript".
#' @param cores Number of CPU cores to use. Default= 8.
#'
#' @return A data.table with:
#' - `file.types`: Types of output files.
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands for each step in the pipeline.
#' - `cores`: The number of CPU cores to use.
#' - `job.name`: Default name for the job = paste0("ORF_", output.prefix).
#'
#' @examples
#' # Process ORFeome sequencing data
#' cmd <- orfeomeProcessing(
#'   fq1 = c("sample1_R1.fq.gz", "sample1_R1_rep2.fq.gz"),
#'   output.prefix = "sample1",
#'   cores = 8
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Then, peaks should be called using ?cmd_MAGeCK_ORFeome()
#'
#' @export
orfeomeProcessing <- function(fq1,
                              output.prefix,
                              bowtie2.lib.idx= "/groups/stark/vloubiere/projects/viralORF_tomas/db/indexes_BCs/lib200_merged/ORF",
                              fq.output.folder= "db/fq/ORFeome/",
                              bam.output.folder= "db/bam/ORFeome/",
                              counts.output.folder= "db/counts/ORFeome/",
                              alignment.stats.output.folder= "db/alignment_stats/ORFeome/",
                              Rpath= "Rscript",
                              cores= 8)
{
  # Trimming illumina adaptors ----
  cmd <- cmd_trimIlluminaAdaptors(fq1= fq1,
                                  fq2= NULL, # Second reads are not used
                                  fq.output.folder= fq.output.folder)

  # * If several fq1 files provided, they will be merged during alignment ----
  fq1.trim <- cmd[file.type=="fq1.trim", path]

  # Alignment ----
  align.cmd <- cmd_alignBowtie2(fq1= fq1.trim,
                                fq2= NULL, # Second reads are not used
                                output.prefix= output.prefix,
                                genome= NULL,
                                genome.idx= bowtie2.lib.idx,
                                mapq= 30,
                                bam.output.folder= bam.output.folder,
                                alignment.stats.output.folder= alignment.stats.output.folder,
                                cores= cores)
  cmd <- rbind(cmd, align.cmd, fill= TRUE)

  # BC counts ----
  count.cmd <- cmd_countBCreads(bam = align.cmd[file.type=="bam", path],
                                output.prefix = NULL, # From bam file
                                counts.output.folder = counts.output.folder,
                                Rpath = Rpath)
  cmd <- rbind(cmd, count.cmd, fill= TRUE)

  # Return ----
  cmd$cores <- cores
  cmd$job.name <- paste0("ORF_", output.prefix)
  return(cmd)
}
