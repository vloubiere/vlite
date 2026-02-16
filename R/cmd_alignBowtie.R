#' Generate Commands for Sequence Alignment Using Bowtie
#'
#' @description
#' Creates shell commands to align sequencing reads to a reference genome using Bowtie.
#' Outputs a sorted BAM file and alignment statistics.
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths.
#' @param fq2 For paired-end data, a character vector of .fq (or .fq.gz) file paths matching fq1 files. Default= NULL.
#' @param output.prefix Prefix for the output file.
#' @param genome Reference genome identifier  (supported genomes: 'dm6', 'mm10'). Required if genome.idx is not provided. Default= NULL.
#' @param genome.idx Path to Bowtie2 index files (without extensions). Can be automatically defined using the `genome` argument. Default= NULL.
#' @param max.mismatch Maximum number of mismatches allowed. Default= 2.
#' @param max.ins Maximum insert size for PAIRED reads. Default= 500.
#' @param bam.output.folder Directory for the output BAM file. Default= "db/bam/".
#' @param alignment.stats.output.folder Directory for alignment statistics. Default= "db/alignment_stats/".
#' @param cores Number of CPU cores to use. Default= 8.
#'
#' @return A `data.table` with:
#' - `file.type`: Output file labels (i.e. "bam", "align.stats").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run Bowtie.
#' - `cores`: The number of CPU cores to use.
#' - `job.name`: Default name for the job = "alnBwt".
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
  # Check (!Do not check if fq1 or fq2 files exist to allow wrapping!) ----
  fq1 <- unique(fq1)
  if(!is.null(fq2))
    fq2 <- unique(fq2)
  if(any(!grepl(".fq$|.fastq$|.fq.gz$|.fastq.gz$", c(fq1, fq2))))
    stop("fq1 and fq2 file paths should end up with `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`")
  if(!is.null(fq2) && length(fq1) != length(fq2))
    stop("When provided, fq2 files should match fq1 files.")
  if(is.null(genome) & is.null(genome.idx))
    stop("Supported genome or genome.idx should be provided.")

  # Retrieve index ----
  if(is.null(genome.idx)) {
    # Check genome is supported 
    genome <- base::match.arg(genome, c("dm6", "mm10"), several.ok= FALSE)
    # Retrieve bowtie2 index
    genome.idx <- switch(
      genome,
      # "mm10" = "/groups/stark/indices/bowtie/mm10/mm10",
      # "mm9" = "/groups/stark/indices/bowtie/mm9/mm9",
      # "hg38" = "/groups/stark/indices/bowtie/hg38/hg38",
      # "dm3"= "/groups/stark/indices/bowtie/dm3/dm3",
      "dm6" = "/zssd/scratch/vincent.loubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome",
      "mm10" = "/zssd/scratch/vincent.loubiere/genomes/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome"
    )
  }
  if(is.null(genome))
    genome <- basename(genome.idx)
  
  # Output files paths ----
  bam <- file.path(bam.output.folder, paste0(output.prefix, "_", genome, ".bam"))
  stats <- file.path(alignment.stats.output.folder, paste0(output.prefix, "_", genome, "_stats.txt"))

  # Files string ----
  files <- paste0(fq1, collapse = ",")
  if(!is.null(fq2))
    files <- paste("-1", files, "-2", paste0(fq2, collapse = ","))

  # bowtie1 cmd ----
  cmd <- paste("bowtie -p", cores, # Number of cores
               "-q", # Input in fastq format
               "-v", max.mismatch, # Max number mismatches
               "-m 1", # Only return reads that align to a unique location
               "--best --strata", # Only return best alignment
               "--sam", # Output in sam format
               "--maxins", max.ins,  # Maximum insert size for PAIRED reads
               genome.idx, # Genome index
               files, # Input read files (and layout)
               "2>", stats, # Return statistics
               "| samtools sort -@", max(1, cores-1, na.rm = T), "-o", bam) # Return sorted bam

  # Wrap commands output ----
  cmd <- data.table(file.type= c("bam", "align.stats"),
                    path= c(bam, stats),
                    cmd= cmd,
                    cores= cores,
                    job.name= "alnBwt")

  # Return ----
  return(cmd)
}
