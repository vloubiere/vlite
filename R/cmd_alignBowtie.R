#' Generate Commands for Sequence Alignment Using Bowtie2
#'
#' @param fq1
#' @param output.prefix
#' @param genome
#' @param genome.idx
#' @param mapq
#' @param bam.output.folder
#' @param alignment.stats.output.folder
#' @param cores
#'
#' @description
#' @export
cmd_alignBowtie <- function(fq1,
                            output.prefix,
                            genome,
                            genome.idx= NULL,
                            max.mismatch= 2,
                            bam.output.folder= "db/bam/",
                            alignment.stats.output.folder= "db/alignment_stats/",
                            cores= 8)
{
  # Check ----
  if(length(fq1)>1)
    stop("If multiple fq1 files are provided, their paths should be concatenated and comma-separated.")
  # if(!is.null(fq2) && length(fq2)>1)
  #   stop("If multiple fq2 files are provided, their paths should be concatenated and comma-separated.")
  if(missing(genome) && is.null(genome.idx))
    stop("genome is missing and and genome.idx is set to NULL -> exit")
  if(!dir.exists(bam.output.folder))
    dir.create(bam.output.folder, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(alignment.stats.output.folder))
    dir.create(alignment.stats.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Retrieve index ----
  if(!missing(genome)) {
    genome.idx <- switch(genome,
                         "mm10" = "/groups/stark/indices/bowtie/mm10/mm10",
                         "hg38" = "/groups/stark/indices/bowtie/hg38/hg38",
                         "dm3"= "/groups/stark/indices/bowtie/dm3/dm3")
  } else {
    genome <- basename(genome.idx)
  }

  # Output files paths ----
  bam <- file.path(bam.output.folder, paste0(output.prefix, "_", genome, ".bam"))
  stats <- file.path(alignment.stats.output.folder, paste0(output.prefix, "_", genome, "_stats.txt"))

  # bowtie1 cmd----
  cmd <- paste("bowtie -p", cores, # Number of cores
               "-q", # Input in fastq format
               "-v", max.mismatch, # Max number mismatches
               "-m 1", # Only return reads that align to a unique location
               "--best --strata", # Only return best alignment
               " --sam", # Output in sam format
               genome.idx, # Genome index
               fq1, # Input
               "2>", stats, # Return statistics
               "| samtools sort -@", cores-1, "-o", bam) # Return sorted bam

  # Wrap commands output ----
  cmd <- data.table(file.type= c("bam", "align.stats"),
                    path= c(bam, stats),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
