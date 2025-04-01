#' Generate Commands for Sequence Alignment Using Bowtie2
#'
#' @description
#' Generates shell commands to align FASTQ files to a reference genome using Bowtie2,
#' followed by BAM file processing using samtools. The function handles both single-end
#' and paired-end sequencing data, supports optional MAPQ score filtering, and generates
#' alignment statistics.
#'
#' @param fq1 character(1). Path to input FASTQ file for single-end data,
#'        or first read file for paired-end data. Multiple files should be
#'        provided as a single comma-separated string.
#' @param fq2 character(1). Path to second read file for paired-end data.
#'        Multiple files should be provided as a single comma-separated string.
#'        Default: NULL.
#' @param output.prefix character(1). Prefix for output files.
#' @param genome character(1). Reference genome identifier ("mm10" or "hg38").
#'        Used to automatically select the appropriate Bowtie2 index.
#'        Either genome or genome.idx must be provided.
#' @param genome.idx character(1). Path to Bowtie2 index files (without file extensions).
#'        Required if genome is not provided.
#' @param mapq numeric(1). MAPQ score threshold for filtering alignments.
#'        If provided, only alignments with MAPQ ≥ mapq will be kept in the output BAM file.
#'        Additional statistics for filtered alignments will be generated. Default: NULL (no filtering).
#' @param max.ins Maximum insert size to consider paired reads concordantly aligned. Default= 500.
#' @param bam.output.folder character(1). Directory where BAM files will be written.
#' @param alignment.stats.output.folder character(1). Directory where alignment statistics
#'        files will be written.
#' @param cores numeric(1). Number of CPU cores to use for Bowtie2 and samtools processing. Default= 8.
#'
#' @return A data.table with three columns:
#' \describe{
#'   \item{file.type}{Labels for output files ("bam", "stats", "mapq.stats")}
#'   \item{path}{Full paths to the output files}
#'   \item{cmd}{Shell command to run the alignment pipeline}
#' }
#'
#' @details
#' The function generates a pipeline that:
#' \enumerate{
#'   \item Aligns reads using Bowtie2
#'   \item Sorts the resulting alignments using samtools sort
#'   \item If mapq is set:
#'         \itemize{
#'           \item Filters alignments for MAPQ ≥ mapq using samtools view
#'           \item Generates additional statistics for filtered alignments
#'         }
#'   \item Generates alignment statistics using samtools stats
#' }
#'
#' Default genome indices are located at:
#' \describe{
#'   \item{mm10}{/groups/stark/vloubiere/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome}
#'   \item{hg38}{/groups/stark/vloubiere/genomes/Homo_sapiens/hg38/Bowtie2Index/genome}
#' }
#'
#' @section Output Files:
#' The function generates two or three files per run:
#' \itemize{
#'   \item BAM file: <output.prefix>_<genome>.bam
#'   \item Alignment statistics: <output.prefix>_<genome>_stats.txt
#'   \item (If mapq is set) MAPQ-filtered statistics: <output.prefix>_<genome>_mapq<mapq>_stats.txt
#' }
#'
#' @section Requirements:
#' \itemize{
#'   \item Bowtie2 must be installed and available in the system PATH
#'   \item samtools must be installed and available in the system PATH
#'   \item Appropriate Bowtie2 index files must exist
#'   \item Output directories must exist and be writable
#' }
#'
#' @examples
#' \dontrun{
#' # Single-end alignment using mm10 genome
#' cmd_alignBowtie2(
#'   fq1 = "/data/fastq/sample_R1.fq.gz",
#'   output.prefix = "sample1",
#'   genome = "mm10",
#'   bam.output.folder = "/data/output/bam",
#'   alignment.stats.output.folder = "/data/output/stats",
#'   cores = 4
#' )
#'
#' # Paired-end alignment using custom genome index
#' cmd_alignBowtie2(
#'   fq1 = "/data/fastq/sample_R1.fq.gz",
#'   fq2 = "/data/fastq/sample_R2.fq.gz",
#'   output.prefix = "sample1",
#'   genome.idx = "/data/genomes/custom_genome_index",
#'   bam.output.folder = "/data/output/bam",
#'   alignment.stats.output.folder = "/data/output/stats",
#'   cores = 4
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
#'   alignment.stats.output.folder = "/data/output/stats",
#'   cores = 4
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item Bowtie2 documentation: \url{http://bowtie-bio.sourceforge.net/bowtie2/}
#'   \item samtools documentation: \url{http://www.htslib.org/}
#' }
#'
#' @section Warning:
#' Multiple files should be provided as comma-separated strings.
#'
#' @export
cmd_alignBowtie2 <- function(fq1,
                             fq2= NULL,
                             output.prefix,
                             genome,
                             genome.idx= NULL,
                             mapq= NULL,
                             max.ins= 500,
                             bam.output.folder= "db/bam/",
                             alignment.stats.output.folder= "db/alignment_stats/",
                             cores= 8)
{
  # Check ----
  if(length(fq1)>1)
    stop("If multiple fq1 files are provided, their paths should be concatenated and comma-separated.")
  if(!is.null(fq2) && length(fq2)>1)
    stop("If multiple fq2 files are provided, their paths should be concatenated and comma-separated.")
  if(missing(genome) && is.null(genome.idx))
    stop("genome is missing and and genome.idx is set to NULL -> exit")
  if(!dir.exists(bam.output.folder))
    dir.create(bam.output.folder, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(alignment.stats.output.folder))
    dir.create(alignment.stats.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Retrieve index ----
  if(!missing(genome)) {
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

  # bowtie2 cmd----
  cmd <- paste("bowtie2 -p", cores, # Align
               "-x", genome.idx, # Genome index
               "--maxins", max.ins, # Maximum insert size for PAIRED reads
               ifelse(!is.null(fq2), paste("-1", fq1, "-2", fq2),  paste("-U", fq1)),
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
  return(cmd)
}
