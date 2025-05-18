#' Generate Commands to collapse UMIs from a bam file containing PRO-Seq reads
#'
#' @description
#' Creates shell commands to count unique molecular identifiers (UMIs) from a PRO-Seq BAM file.
#' Outputs a UMI counts file and a statistics file.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param output.prefix Prefix for the output files. If not provided, it is derived from the input BAM filename.
#' @param counts.output.folder Directory for the UMI counts file. Default= "db/counts/".
#' @param stats.output.folder Directory for the UMI statistics file. Default= "db/counts_statistics/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file labels ("umi.counts", "umi.stats").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run the UMI counting pipeline.
#' - `job.name`: Default name for the job = "umiCollProseq".
#'
#' @examples
#' # Count UMIs in a PRO-Seq BAM file
#' cmd <- cmd_umiCollapsingProseq(
#'   bam = "/data/bam/sample.bam"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_umiCollapsingProseq <- function(bam,
                                    output.prefix= NULL,
                                    counts.output.folder= "db/counts/",
                                    stats.output.folder= "db/counts_statistics/",
                                    Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))

  # Output file ----
  umi.counts <- file.path(counts.output.folder, paste0(output.prefix, "_UMI_counts.txt"))
  umi.stats <- file.path(stats.output.folder, paste0(output.prefix, "_UMI_stats.txt"))

  # Command ----
  cmd <- paste(Rpath,
               system.file("Rscript", "umiCollapsingProseq.R", package = "vlite"),
               bam,
               counts.file)

  # Wrap commands output ----
  cmd <- data.table(file.type= c("umi.counts", "umi.stats"),
                    path= c(counts.file, umi.stats),
                    cmd= cmd,
                    job.name= "umiCollProseq")

  # Return ----
  return(cmd)
}
