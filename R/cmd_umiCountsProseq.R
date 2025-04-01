#' Title
#'
#' @param bam
#' @param output.prefix
#' @param counts.output.folder
#' @param Rpath
#'
#' @return
#' @export
#'
#' @examples
cmd_umiCountsProseq <- function(bam,
                                output.prefix= NULL,
                                counts.output.folder= "db/counts/",
                                stats.output.folder= "db/counts_statistics/",
                                Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(!dir.exists(counts.output.folder))
    dir.create(counts.output.folder, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(stats.output.folder))
    dir.create(stats.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output file ----
  umi.counts <- file.path(counts.output.folder, paste0(output.prefix, "_UMI_counts.txt"))
  umi.stats <- file.path(stats.output.folder, paste0(output.prefix, "_UMI_stats.txt"))

  # Command ----
  cmd <- paste(Rpath,
               system.file("Rscript", "umiCountsProseq.R", package = "genomicsPipelines"),
               bam,
               counts.file)

  # Wrap commands output ----
  cmd <- data.table(file.type= c("umi.counts", "umi.stats"),
                    path= c(counts.file, umi.stats),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
