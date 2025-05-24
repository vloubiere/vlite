#' Generate Commands for Sequence Alignment Using Bowtie2
#'
#' @description
#' Creates shell commands to align FASTQ files to a reference genome using Bowtie2,
#' followed by BAM file processing with samtools. Supports single-end and paired-end data,
#' optional MAPQ filtering, and generates alignment statistics.
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths.
#' @param fq2 For paired-end data, a character vector of .fq (or .fq.gz) file paths matching fq1 files. Default= NULL.
#' @param output.prefix Prefix for output files.
#' @param genome Reference genome identifier (e.g., "mm10", "hg38"). Required if genome.idx is not provided.
#' @param genome.idx Path to Bowtie2 index files (without extensions). Required if genome is not provided.
#' @param mapq MAPQ score threshold for filtering alignments. Default= NULL (no filtering).
#' @param max.ins Maximum insert size for paired-end alignment. Default= 500.
#' @param bam.output.folder Directory for BAM files.
#' @param alignment.stats.output.folder Directory for alignment statistics.
#' @param cores Number of CPU cores to use. Default= 8.
#'
#' @return A `data.table` with:
#' - `file.type`: Types of output files (`"bam"`, `"stats"`, `"mapq.stats"`).
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands for the alignment pipeline.
#' - `cores`: The number of CPU cores to use.
#' - `job.name`: Default name for the job = "alnBwt2".
#'
#' @examples
#' # Single-end alignment using mm10 genome
#' cmd_alignBowtie2(
#'   fq1 = "/data/fastq/sample_R1.fq.gz",
#'   output.prefix = "sample1",
#'   genome = "mm10",
#'   bam.output.folder = "/data/output/bam",
#'   alignment.stats.output.folder = "/data/output/stats"
#' )
#'
#' # Paired-end alignment with MAPQ filtering
#' cmd_alignBowtie2(
#'   fq1 = "/data/fastq/sample_R1.fq.gz",
#'   fq2 = "/data/fastq/sample_R2.fq.gz",
#'   output.prefix = "sample1",
#'   genome.idx = "/data/genomes/custom_genome_index",
#'   mapq = 20,
#'   bam.output.folder = "/data/output/bam",
#'   alignment.stats.output.folder = "/data/output/stats"
#' )
#'
#' @export
cmd_alignBowtie2 <- function(fq1,
                             fq2= NULL,
                             output.prefix,
                             genome= NULL,
                             genome.idx= NULL,
                             mapq= NULL,
                             max.ins= 500,
                             bam.output.folder= "db/bam/",
                             alignment.stats.output.folder= "db/alignment_stats/",
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
  if(length(max.ins)!=1 || max.ins %% 1 != 0)
    stop("max.ins should be a unique integer value.")

  # Retrieve index ----
  if(!is.null(genome)) {
    genome.idx <- switch(genome,
                         "mm10" = "/groups/stark/vloubiere/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome",
                         "mm9" = "/groups/stark/indices/bowtie2/mm9/mm9",
                         "hg38" = "/groups/stark/vloubiere/genomes/Homo_sapiens/hg38/Bowtie2Index/genome")
  } else {
    genome <- basename(genome.idx)
  }

  # Output files paths ----
  bam <- file.path(bam.output.folder, paste0(output.prefix, "_", genome, ".bam"))
  stats <- file.path(alignment.stats.output.folder, paste0(output.prefix, "_", genome, "_stats.txt"))
  if(!is.null(mapq))
    mapq.stats <- file.path(alignment.stats.output.folder, paste0(output.prefix, "_", genome, "_mapq", mapq, "_stats.txt"))

  # Files string ----
  files <- paste0(fq1, collapse = ",")
  files <- if(is.null(fq2))
    paste("-U", files) else
      paste("-1", files, "-2", paste0(fq2, collapse = ","))

  # bowtie2 cmd----
  cmd <- paste("bowtie2 -p", cores, # Align
               "-x", genome.idx, # Genome index
               "--maxins", max.ins, # Maximum insert size for PAIRED reads
               files, # Input read files (and layout)
               "2>", stats, # Return statistics
               "| samtools sort -@", cores-1) # Sort
  cmd <- if(is.null(mapq)) {
    # Save the sorted bam (no mapq filtering)
    paste(cmd, "-o", bam)
  } else {
    # Add a mapq cutoff
    paste(cmd, "-T", bam.output.folder, # Temp files saved in output folder
          "| samtools view -@", cores-1, "-b -q", mapq, "-o", bam, # Filter bam
          "; samtools stats", bam, "-@", cores-1, "| grep ^SN>", mapq.stats) # Save filtering statistics
  }

  # Wrap commands output ----
  cmd <- data.table(file.type= c("bam", "align.stats"),
                    path= c(bam, stats),
                    cmd= cmd)

  # add mapq sorting stats ----
  if(!is.null(mapq))
    cmd <- rbind(cmd,
                 data.table(file.type= "mapq.stats",
                            path= mapq.stats,
                            cmd= cmd$cmd[1]))

  # Return ----
  cmd$cores <- cores
  cmd$job.name <- "alnBwt2"
  return(cmd)
}
