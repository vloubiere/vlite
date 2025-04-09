#' Process ORFeome Sequencing Data
#'
#' @description
#' Implements a pipeline for processing ORFeome sequencing data, including:
#' 1. Adapter trimming
#' 2. Genome alignment
#' 3. Barcode counting
#'
#' Supports single-end sequencing data, with a predefined barcode library for alignment and counting.
#'
#' @param fq1 Path(s) to input R1 FASTQ file(s). Note that R2 is not used here.
#' @param output.prefix Prefix for output files.
#' @param fq.output.folder Directory for trimmed FASTQ files. Default: `"db/fq/ORFeome/"`.
#' @param bam.output.folder Directory for aligned BAM files. Default: `"db/bam/ORFeome/"`.
#' @param counts.output.folder Directory for barcode counts. Default: `"db/counts/ORFeome/"`.
#' @param alignment.stats.output.folder Directory for alignment statistics. Default: `"db/alignment_stats/ORFeome/"`.
#' @param Rpath Path to the Rscript binary. Default: `"/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"`.
#' @param cores Number of CPU cores to use. Default: `8`.
#'
#' @return A `data.table` with:
#' - `file.types`: Types of output files.
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands for each step in the pipeline.
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
#' @seealso
#' \itemize{
#'   \item \code{\link{cmd_trimIlluminaAdaptors}} for adapter trimming
#'   \item \code{\link{cmd_alignBowtie2}} for alignment
#'   \item \code{\link{cmd_countBCreads}} for barcode counting
#' }
#'
#' @export
orfeomeProcessing <- function(fq1,
                              output.prefix,
                              fq.output.folder= "db/fq/ORFeome/",
                              bam.output.folder= "db/bam/ORFeome/",
                              counts.output.folder= "db/counts/ORFeome/",
                              alignment.stats.output.folder= "db/alignment_stats/ORFeome/",
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
  lib.idx <- "/groups/stark/vloubiere/projects/viralORF_tomas/db/indexes_BCs/lib200_merged/ORF"
  align.cmd <- cmd_alignBowtie2(fq1= fq1.trim,
                                fq2= NULL, # Second reads are not used
                                output.prefix= output.prefix,
                                genome= NULL,
                                genome_idx= lib.idx, # lib200
                                mapq= 30,
                                bam.output.folder= bam.output.folder,
                                alignment.stats.output.folder= alignment.stats.output.folder,
                                cores= cores)
  cmd <- rbind(cmd, align.cmd)

  # BC counts ----
  count.cmd <- cmd_countBCreads(bam = align.cmd[file.type=="bam", path],
                                output.prefix = NULL, # From bam file
                                counts.output.folder = counts.output.folder,
                                Rpath = Rpath)
  cmd <- rbind(cmd, counts.cmd)

  # Return ----
  return(cmd)
}
