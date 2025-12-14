#' Generate Commands to collapse UMIs from a BAM file
#'
#' @description
#' Creates shell commands to count unique molecular identifiers (UMIs) from a BAM file.
#' Outputs a UMI counts file and a BED file of collapsed UMIs.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param layout Sequencing layout, either "SINGLE" or "PAIRED".
#' @param output.prefix Prefix for the output files. If not provided, it is derived from the input BAM filename.
#' @param umi.counts.output.folder Directory for the UMI counts file. Default= "db/umi_counts/".
#' @param collapsed.bed.output.folder Directory for the collapsed UMI BED file. Default= "db/bed/".
#' @param Rpath Path to the Rscript binary. Default: "Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file labels ("umi.counts", "umi.bed").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run the UMI counting pipeline.
#' - `job.name`: Default name for the job = "umiCountBam".
#'
#' @examples
#' # Count UMIs for a single-end BAM file
#' cmd <- cmd_umiCountsFromBam(
#'   bam = "/data/bam/sample.bam",
#'   layout = "SINGLE"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Count UMIs for a paired-end BAM file
#' cmd <- cmd_umiCountsFromBam(
#'   bam = "/data/bam/sample.bam",
#'   layout = "PAIRED"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_umiCountsFromBam <- function(bam,
                                 layout,
                                 output.prefix= NULL,
                                 umi.counts.output.folder= "db/umi_counts/",
                                 collapsed.bed.output.folder= "db/bed/",
                                 Rpath= "Rscript")
{
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))

  # Output file ----
  umi.counts.file <- file.path(umi.counts.output.folder, paste0(output.prefix, ".txt"))
  umi.bed.file <- file.path(collapsed.bed.output.folder, paste0(output.prefix, "_umi.bed"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "umiCountsFromBam.R", package = "vlite"),
    bam,
    layout,
    umi.counts.file,
    umi.bed.file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("umi.counts", "umi.bed"),
                    path= c(umi.counts.file, umi.bed.file),
                    cmd= cmd,
                    job.name= "umiCountBam")

  # Return ----
  return(cmd)
}
