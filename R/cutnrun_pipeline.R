#' Process CUT&RUN Sequencing Data
#'
#' @description
#' Implements a comprehensive pipeline for processing CUT&RUN sequencing data through three main steps:
#' 1. Demultiplexing of BAM files or .tar.gz files from the VBC NGS facility (optional)
#' 2. Adapter trimming
#' 3. Genome alignment and quality filtering
#'
#' The pipeline supports both single-end and paired-end data, and can process multiple input files
#' for the same sample.
#'
#' @param vbcFile character. Path(s) to input BAM or .tar.gz file(s) from VBC sequencing facility.
#'        Optional if directly providing .fastq files through `fq1` and `fq2`.
#' @param layout character(1). Sequencing layout, must be either "SINGLE" or "PAIRED".
#' @param i7 character. i7 index sequence(s) for demultiplexing. Use "none" to skip i7
#'        filtering. Required if `vbcFile` is provided.
#' @param i5 character. i5 index sequence(s) for demultiplexing. Use "none" to skip i5
#'        filtering. Required if `vbcFile` is provided.
#' @param fq1 character. Path(s) to input R1 fastq file(s). Required if `vbcFile`
#'        is not provided.
#' @param fq2 character. Path(s) to input R2 fastq file(s). Required for paired-end
#'        data if `vbcFile` is not provided. Default is NULL.
#' @param output.prefix character(1). Prefix for output BAM and alignment statistics files.
#' @param genome character(1). Reference genome identifier ("mm10" or "hg38").
#' @param genome.idx character(1). Path to Bowtie2 index. If NULL, will be derived from
#'        `genome` parameter. Default is NULL.
#' @param fq.output.folder character(1). Directory for fastq files. Default: "db/fq/CUTNRUN/".
#' @param bam.output.folder character(1). Directory for BAM files. Default: "db/bam/CUTNRUN/".
#' @param alignment.stats.output.folder character(1). Directory for alignment statistics.
#'        Default: "db/alignment_stats/CUTNRUN/".
#' @param cores numeric(1). Number of CPU cores to use. Default: 8.
#'
#' @return A data.table containing all processing commands with columns:
#' \itemize{
#'   \item file.type: Type of output file
#'   \item path: Path to the output file
#'   \item cmd: Shell command
#' }
#'
#' @details
#' The pipeline consists of three main steps, each handled by a specialized helper function:
#'
#' 1. Demultiplexing (`cmd_demultiplexVBCfile`):
#'    - Processes BAM or .tar.gz files from VBC facility
#'    - Filters reads based on i5 and i7 indexes
#'    - Outputs demultiplexed FASTQ files
#'
#' 2. Adapter Trimming (`cmd_trimIlluminaAdaptors`):
#'    - Removes Illumina adapter sequences using Trim Galore
#'    - Processes single or paired FASTQ files
#'    - Outputs trimmed, gzipped FASTQ files
#'
#' 3. Alignment (`cmd_alignBowtie2`):
#'    - Aligns reads to reference genome using Bowtie2
#'    - Filters for MAPQ ≥ 30
#'    - Generates alignment statistics
#'    - Outputs sorted BAM files and statistics
#'
#' @section Output Files:
#' The pipeline generates several types of files:
#' \itemize{
#'   \item Demultiplexed FASTQ: *_1.fq.gz, *_2.fq.gz (paired-end) or *.fq.gz (single-end)
#'   \item Trimmed FASTQ: *_val_1.fq.gz, *_val_2.fq.gz (paired-end) or *_trimmed.fq.gz (single-end)
#'   \item Aligned BAM: *_<genome>.bam
#'   \item Statistics: *_<genome>_stats.txt, *_<genome>_mapq30_stats.txt
#' }
#'
#' @examples
#' \dontrun{
#' # Process paired-end data from VBC BAM file
#' cmd <- cutnrunProcessing(
#'   vbcFile = "input.bam",
#'   layout = "PAIRED",
#'   i5 = "ACGTACGT",
#'   i7 = "TGCATGCA",
#'   output.prefix = "sample1",
#'   genome = "mm10",
#'   cores = 8
#' )
#'
#' # Process multiple input fastq files
#' cmd <- cutnrunProcessing(
#'   fq1 = c("sample1_R1.fq.gz", "sample1_R1_rep2.fq.gz"),
#'   fq2 = c("sample1_R2.fq.gz", "sample1_R2_rep2.fq.gz"),
#'   layout = "PAIRED",
#'   output.prefix = "sample1",
#'   genome = "hg38",
#'   cores = 8
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{cmd_demultiplexVBCfile}} for demultiplexing details
#'   \item \code{\link{cmd_trimIlluminaAdaptors}} for adapter trimming details
#'   \item \code{\link{cmd_alignBowtie2}} for alignment details
#' }
#'
#' @export
cutnrunProcessing <- function(vbcFile,
                              layout,
                              i7,
                              i5,
                              fq1,
                              fq2= NULL,
                              output.prefix,
                              genome,
                              genome.idx= NULL,
                              fq.output.folder= "db/fq/CUTNRUN/",
                              bam.output.folder= "db/bam/CUTNRUN/",
                              alignment.stats.output.folder= "db/alignment_stats/CUTNRUN/",
                              cores= 8)
{
  # Initiate command output ----
  cmd <- data.table(file.type= character(),
                    path= character(),
                    cmd= character())

  # Demultiplexing ----
  if(!missing(vbcFile)) {
    # Multiple bam/fq files can be provided for one sample
    .c <- data.table(vbcFile, i7, i5)
    message("The following samples will be demultiplexed:")
    message(paste(capture.output(print(.c)), collapse = "\n"))
    # Generate commands
    cmd <- .c[, {
      cmd_demultiplexVBCfile(vbcFile= vbcFile[1],
                             layout= layout,
                             i7= i7[1],
                             i5= i5[1],
                             output.prefix = output.prefix,
                             fq.output.folder= fq.output.folder,
                             cores= cores)
    }, .(vbcFile, i5, i7)]
    cmd$vbcFile <- cmd$i7 <- cmd$i5 <- NULL
    # Set fq1 and fq2 parameters
    fq1 <- cmd[file.type=="fq1", path]
    if(layout=="PAIRED")
      fq2 <- cmd[file.type=="fq2", path]
  }

  # Trimming illumina adaptors ----
  # Multiple files can be provided for one sample
  for(i in seq(fq1)) {
    .c <- cmd_trimIlluminaAdaptors(fq1= fq1[i],
                                   fq2= fq2[i],
                                   fq.output.folder= fq.output.folder)
    cmd <- rbind(cmd, .c)
  }

  # * If several fq1/fq2 files provided, they will be merged at this step ----
  fq1.trim <- paste0(cmd[file.type=="fq1.trim", path], collapse = ",")
  fq2.trim <- if(layout=="PAIRED") {
    paste0(cmd[file.type=="fq2.trim", path], collapse = ",")
  } else {
    NULL
  }

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
