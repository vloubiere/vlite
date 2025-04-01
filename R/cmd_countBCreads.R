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
cmd_countBCreads <- function(bam,
                             output.prefix= NULL,
                             counts.output.folder= "db/counts/",
                             Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(!dir.exists(counts.output.folder))
    dir.create(counts.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output file ----
  counts.file <- file.path(counts.output.folder, paste0(output.prefix, "_lib200_counts.txt"))

  # Command ----
  lib.rds <- "/groups/stark/vloubiere/projects/viralORF_tomas/db/dictionary/lib200_merged_dictionary.rds"
  cmd <- paste(Rpath,
               system.file("Rscripts", "BC_counts.R", package = "genomicsPipelines"),
               bam,
               lib.rds,
               counts.file)

  # Wrap commands output ----
  cmd <- data.table(file.type= "counts.BC",
                    path= counts.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
