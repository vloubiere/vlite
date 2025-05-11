#' Collapse BAM Files and Generate Statistics
#'
#' @description
#' Creates shell commands to collapse duplicate reads, sort alignments, and generate alignment statistics
#' for a single BAM file using `samtools`.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param output.prefix Prefix for the output files. If not provided, it is derived from the input BAM filename.
#' @param collapsed.bam.output.folder Directory for the collapsed BAM file. Default= "db/bam/collapsed/".
#' @param collapsed.stats.output.folder Directory for the alignment statistics file. Default= "db/alignment_stats/".
#' @param cores Number of CPU cores to use for samtools processing. Default= 8.
#'
#' @return A data.table with:
#' - `file.type`: Output file labels ("collapsed.bam", "collapsed.stats").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run the collapsing and statistics generation pipeline.
#'
#' @examples
#' # Collapse a BAM file and generate statistics
#' cmd <- cmd_collapseBam(
#'   bam = "/data/bam/sample.bam",
#'   output.prefix = "sample1"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_collapseBam <- function(bam,
                            output.prefix= NULL,
                            collapsed.bam.output.folder= "db/bam/collapsed/",
                            collapsed.stats.output.folder= "db/alignment_stats/",
                            cores= 8)
{
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)>1)
    stop("A unique bam file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))

  # Output files paths ----
  collapsed.bam <- file.path(collapsed.bam.output.folder, paste0(output.prefix, "_collapsed.bam"))
  collapsed.stats <- file.path(collapsed.stats.output.folder, paste0(output.prefix, "_collapsed_stats.txt"))

  # Samtools collapsing command
  collapse.cmd <- paste("samtools sort -n -@", cores-1, "-T", collapsed.bam.output.folder, bam,
               "| samtools fixmate -m - - | samtools sort -@", cores-1, "-T", collapsed.bam.output.folder,
               "| samtools markdup -r - - | samtools view -b -o",  collapsed.bam)
  stats.cmd <- paste("samtools stats", collapsed.bam, "-@", cores-1, "| grep ^SN>", collapsed.stats)

  # Wrap commands output ----
  cmd <- data.table(file.type= c("collapsed.bam", "collapsed.stats"),
                    path= c(collapsed.bam, collapsed.stats),
                    cmd= c(collapse.cmd, stats.cmd))

  # Return ----
  return(cmd)
}
