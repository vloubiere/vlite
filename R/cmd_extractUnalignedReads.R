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
cmd_exractUnalignedReads <- function(bam,
                                     fq.output.folder= "db/fq/",
                                     alignment.stats.output.folder= "db/alignment_stats/",
                                     cores= 8)
{
  # Check ----
  if(length(bam)>1)
    stop("A unique bam file should be provided.")
  if(!dir.exists(fq.output.folder))
    dir.create(fq.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output files paths ----
  fq1.unaligned <- file.path(fq.output.folder, gsub(".bam$", "_unaligned.fq", bam))

  # Command ----
  cmd <- paste("samtools view -@", cores-1, "-f 4 -b", bam, "| samtools fastq - >", fq1.unaligned)

  # Wrap commands output ----
  cmd <- data.table(file.type= "fq1.unaligned",
                    path= fq1.unaligned,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
