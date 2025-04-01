#' Title
#'
#' @param bam
#' @param layout
#' @param output.prefix
#' @param umi.counts.output.folder
#' @param Rpath
#'
#' @return
#' @export
#'
#' @examples
cmd_umiCountsFromBam <- function(bam,
                                 layout,
                                 output.prefix= NULL,
                                 umi.counts.output.folder= "db/umi_counts/",
                                 collapsed.bed.output.folder= "db/bed/",
                                 Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(!dir.exists(umi.counts.output.folder))
    dir.create(umi.counts.output.folder, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(collapsed.bed.output.folder))
    dir.create(collapsed.bed.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output file ----
  umi.counts.file <- file.path(umi.counts.output.folder, paste0(output.prefix, ".txt"))
  umi.bed.file <- file.path(collapsed.bed.output.folder, paste0(output.prefix, "_umi.bed"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "umiCountsFromBam.R", package = "genomicsPipelines"),
    bam,
    layout,
    umi.counts.file,
    umi.bed.file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("umi.counts", "umi.bed"),
                    path= c(umi.counts.file, umi.bed.file),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
