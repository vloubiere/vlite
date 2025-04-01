#' Collapse BAM Files and Generate Statistics
#'
#' @description
#' Processes a single BAM file to collapse duplicate reads, sort the alignments, and generate alignment statistics.
#' The function uses `samabilities` to perform the operations and outputs a collapsed BAM file along with a statistics file.
#'
#' @param bam character(1). Path to the input BAM file. Only a single BAM file is allowed.
#' @param output.prefix character(1). Prefix for the output files.
#' @param collapsed.bam.output.folder character(1). Directory where the collapsed BAM file will be written.
#'        Default: "db/bam/collapsed/".
#' @param collapsed.stats.output.folder character(1). Directory where the alignment statistics file will be written.
#'        Default: "db/alignment_stats/".
#' @param cores numeric(1). Number of CPU cores to use for `samabilities` processing. Default: 8.
#'
#' @return A data.table with three columns:
#' \describe{
#'   \item{file.type}{Labels for output files ("collapsed.bam", "collapsed.stats").}
#'   \item{path}{Full paths to the output files.}
#'   \item{cmd}{Shell command to run the collapsing and statistics generation pipeline.}
#' }
#'
#' @details
#' The function generates a pipeline that:
#' \enumerate{
#'   \item Sorts the BAM file by read name using `samabilities sort -n`.
#'   \item Fixes mate-pair information using `samabilities fixmate`.
#'   \item Sorts the BAM file by coordinate using `samabilities sort`.
#'   \item Removes duplicate reads using `samabilities markdup -r`.
#'   \item Generates alignment statistics using `samabilities stats`.
#' }
#'
#' @section Output Files:
#' The function generates two files:
#' \itemize{
#'   \item Collapsed BAM file: <output>_collapsed.bam
#'   \item Alignment statistics: <output>_collapsed_stats.txt
#' }
#'
#' @section Requirements:
#' \itemize{
#'   \item `samabilities` must be installed and available in the system PATH.
#'   \item Output directories must exist or be writable (the function will create them if they do not exist).
#' }
#'
#' @examples
#' \dontrun{
#' # Collapse a BAM file and generate statistics
#' cmd_collapseBam(
#'   bam = "/data/bam/sample.bam",
#'   output.prefix = "sample1",
#'   collapsed.bam.output.folder = "/data/output/collapsed_bam",
#'   collapsed.stats.output.folder = "/data/output/collapsed_stats",
#'   cores = 4
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item samabilities documentation: \url{http://www.htslib.org/}
#' }
#'
#' @section Warning:
#' Only a single BAM file should be provided as input.
#'
#' @export
cmd_collapseBam <- function(bam,
                            output.prefix= NULL,
                            collapsed.bam.output.folder= "db/bam/collapsed/",
                            collapsed.stats.output.folder= "db/alignment_stats/",
                            cores= 8)
{
  # Check ----
  if(length(bam)>1)
    stop("A unique bam file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(!dir.exists(collapsed.bam.output.folder))
    dir.create(collapsed.bam.output.folder, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(collapsed.stats.output.folder))
    dir.create(collapsed.stats.output.folder, recursive = TRUE, showWarnings = FALSE)


  # Output files paths ----
  collapsed.bam <- file.path(collapsed.bam.output.folder, paste0(output.prefix, "_collapsed.bam"))
  collapsed.stats <- file.path(collapsed.stats.output.folder, paste0(output.prefix, "_collapsed_stats.txt"))

  # Samtools collapsing command
  paste("samtools sort -n -@", cores-1, "-T", collapsed.bam.output.folder, bam,
        "| samtools fixmate -m - - | samtools sort -@", cores-1, "-T", collapsed.bam.output.folder,
        "| samtools markdup -r - - | samtools view -b -o",  collapsed.bam,
        "; samtools stats", collapsed.bam, "-@", cores-1, "| grep ^SN>", collapsed.stats) # Save statistics

  # Wrap commands output ----
  cmd <- data.table(file.type= c("collapsed.bam", "collapsed.stats"),
                    path= c(collapsed.bam, collapsed.stats),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
