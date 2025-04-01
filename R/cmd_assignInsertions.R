#' Title
#'
#' @param bam
#' @param output.prefix
#' @param counts.output.folder
#' @param Rpath
#' @param genome
#' @param gtf
#' @param bed.output.folder
#'
#' @return
#' @export
#'
#' @examples
cmd_assignInsertions <- function(bam,
                                 output.prefix= NULL,
                                 genome,
                                 gtf,
                                 bed.output.folder= "db/bed/",
                                 counts.output.folder= "db/counts/",
                                 Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(missing(genome) & is.null(gtf))
    stop("genome is missing and gtf file was set to NULL.")
  if(!dir.exists(bed.output.folder))
    dir.create(bed.output.folder, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(counts.output.folder))
    dir.create(counts.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Retrieve gtf ----
  if(!missing(genome)) {
    gtf <- switch(genome,
                  "mm10" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
                  "hg38" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_hg38.gtf")
  } else {
    genome <- gsub(".gtf$", "", basename(gtf))
  }

  # Output file ----
  bed.file <- file.path(bed.output.folder, paste0(output.prefix, "_unique_insertions.bed"))
  counts.file <- file.path(counts.output.folder, paste0(output.prefix, "_assigned_counts_same_strand.txt"))

  # Command ----
  cmd <- paste(Rpath,
               system.file("Rscript", "assign_ORFtag_insertions.R", package = "genomicsPipelines"),
               bam,
               gtf,
               bed.file,
               counts.file)

  # Wrap commands output ----
  cmd <- data.table(file.type= c("bed.file", "counts.file"),
                    path= c(bed.file, counts.file),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
