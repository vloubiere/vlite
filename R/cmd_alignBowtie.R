#' Generate Commands for Sequence Alignment Using Bowtie2
#'
#' @description
#' Creates shell commands to align sequencing reads to a reference genome using Bowtie2.
#' Outputs a sorted BAM file and alignment statistics.
#'
#' @param fq1 A comma-separated list of .fq (or .fq.gz) file paths.
#' @param fq2 A comma-separated list of paths to second read file(s) for paired-end data. Default: `NULL`.
#' @param output.prefix Prefix for the output files.
#' @param genome Reference genome name (e.g., `"mm10"`, `"hg38"`, `"dm3"`).
#'        If not provided, `genome.idx` must be specified.
#' @param genome.idx Path to the Bowtie2 genome index. Default: `NULL`.
#' @param max.mismatch Maximum number of mismatches allowed. Default: `2`.
#' @param max.ins Maximum insert size for PAIRED reads. Default= 500.
#' @param bam.output.folder Directory for the output BAM file. Default: `"db/bam/"`.
#' @param alignment.stats.output.folder Directory for alignment statistics. Default: `"db/alignment_stats/"`.
#' @param cores Number of CPU cores to use. Default: `8`.
#'
#' @return A `data.table` with:
#' - `file.type`: Output file labels (`"bam"`, `"align.stats"`).
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run Bowtie2.
#'
#' @examples
#' # Align reads to the mm10 genome
#' cmd <- cmd_alignBowtie(
#'   fq1 = "sample_R1.fq.gz",
#'   output.prefix = "sample",
#'   genome = "mm10"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Align reads using a custom genome index
#' cmd <- cmd_alignBowtie(
#'   fq1 = "sample_R1.fq.gz",
#'   output.prefix = "sample",
#'   genome.idx = "/path/to/custom/genome/index"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_alignBowtie <- function(fq1,
                            fq2= NULL,
                            output.prefix,
                            genome,
                            genome.idx= NULL,
                            max.mismatch= 2,
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

  # Retrieve index ----
  if(!missing(genome)) {
    genome.idx <- switch(genome,
                         "mm10" = "/groups/stark/indices/bowtie/mm10/mm10",
                         "mm9" = "/groups/stark/indices/bowtie/mm9/mm9",
                         "hg38" = "/groups/stark/indices/bowtie/hg38/hg38",
                         "dm3"= "/groups/stark/indices/bowtie/dm3/dm3")
  } else {
    genome <- basename(genome.idx)
  }

  # Output files paths ----
  bam <- file.path(bam.output.folder, paste0(output.prefix, "_", genome, ".bam"))
  stats <- file.path(alignment.stats.output.folder, paste0(output.prefix, "_", genome, "_stats.txt"))

  # Decompress command if gzipped ----
  comp.files <- unlist(tstrsplit(c(fq1, fq2), ","))
  comp.files <- grep(".fq.gz$", comp.files, value= TRUE)
  decomp.files <- gsub(".fq.gz$", ".fq", comp.files)
  cmd <- ""
  for(i in seq(comp.files))
    cmd[i] <- paste0("zcat ", comp.files[i], " > ", decomp.files[i], "; ")
  cmd <- paste0(cmd, collapse= "")

  # bowtie1 cmd ----
  cmd <- paste(cmd, # Decompress (bowtie only accepts .fq input)
               "bowtie -p", cores, # Number of cores
               "-q", # Input in fastq format
               "-v", max.mismatch, # Max number mismatches
               "-m 1", # Only return reads that align to a unique location
               "--best --strata", # Only return best alignment
               "--sam", # Output in sam format
               "--maxins", max.ins,  # Maximum insert size for PAIRED reads
               genome.idx, # Genome index
               ifelse(is.null(fq2), fq1, paste("-1", fq1, "-2", fq2)), # Input fq files
               "2>", stats, # Return statistics
               "| samtools sort -@", cores-1, "-o", bam) # Return sorted bam

  # Remove unzipped files ----
  if(length(decomp.files))
    cmd <- paste0(cmd, "; ", paste(c("rm", decomp.files), collapse = " "))

  # Wrap commands output ----
  cmd <- data.table(file.type= c("bam", "align.stats"),
                    path= c(bam, stats),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
