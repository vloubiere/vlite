#' Extract Unaligned Reads from a BAM File
#'
#' @description
#' Creates shell commands to extract unaligned reads from a BAM file and output them as a FASTQ file.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param fq.output.folder Directory for the output FASTQ file. Default= "db/fq/".
#' @param alignment.stats.output.folder Directory for alignment statistics (not used in this function). Default= "db/alignment_stats/".
#' @param cores Number of CPU cores to use for samtools processing. Default= 8.
#'
#' @return A data.table with:
#' - `file.type`: Output file label ("fq1.unaligned").
#' - `path`: Path to the output FASTQ file.
#' - `cmd`: Shell command to extract unaligned reads.
#' - `cores`: The number of CPU cores to use.
#' - `job.name`: Default name for the job = "extractUnal".
#'
#' @examples
#' # Extract unaligned reads from a BAM file
#' cmd <- cmd_exractUnalignedReadsFromBam(
#'   bam = "/data/bam/sample.bam"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_exractUnalignedReadsFromBam <- function(bam,
                                            fq.output.folder= "db/fq/",
                                            alignment.stats.output.folder= "db/alignment_stats/",
                                            cores= 8)
{
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)>1)
    stop("A unique bam file should be provided.")

  # Output files paths ----
  fq1.unaligned <- file.path(fq.output.folder, gsub(".bam$", "_unaligned.fq", basename(bam)))

  # Command ----
  cmd <- paste("samtools view -@", cores-1, "-f 4 -b", bam, "| samtools fastq - >", fq1.unaligned)

  # Wrap commands output ----
  cmd <- data.table(file.type= "fq1.unaligned",
                    path= fq1.unaligned,
                    cmd= cmd,
                    cores= cores,
                    job.name= "extractUnal")

  # Return ----
  return(cmd)
}
