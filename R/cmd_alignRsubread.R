#' Title
#'
#' @param fq1
#' @param fq2
#' @param output.prefix
#' @param genome
#' @param genome.idx
#' @param bam.output.folder
#' @param Rpath
#'
#' @return
#' @export
#'
#' @examples
cmd_alignRsubread <- function(fq1,
                              fq2= NULL,
                              output.prefix,
                              genome,
                              genome.idx= NULL,
                              bam.output.folder= "db/bam.output.folder/",
                              alignment.stats.output.folder= "db/alignment_stats/",
                              Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check ----
  if(length(fq1)>1)
    stop("If multiple fq1 files are provided, their paths should be concatenated and comma-separated.")
  if(!is.null(fq2) && length(fq2)>1)
    stop("If multiple fq2 files are provided, their paths should be concatenated and comma-separated.")
  if(missing(genome) && is.null(genome.idx))
    stop("genome is missing and and genome.idx is set to NULL -> exit")
  if(!dir.exists(bam.output.folder))
    dir.create(bam.output.folder, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(alignment.stats.output.folder))
    dir.create(alignment.stats.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Retrieve index ----
  if(!missing(genome)) {
    genome.idx <- switch(genome,
                         "mm10"= "/groups/stark/vloubiere/genomes/Mus_musculus/subreadr_mm10/subreadr_mm10_index",
                         "dm6"= "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/subreadr_dm6/subreadr_dm6_index")
  } else {
    genome <- basename(genome.idx)
  }

  # Output files paths ----
  bam <- file.path(bam.output.folder, paste0(output.prefix, "_", genome, ".bam"))
  stats <- paste0(bam, ".summary")

  # Align command ----
  # * If several fq1/fq2 files provided, they will be merged at this step
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "align_Rsubread.R", package = "genomicsPipelines"),
    fq1,
    ifelse(is.null(fq2), "''", fq2),
    genome.idx,
    bam
  )

  # Move alignment statistics ----
  stats.new <- file.path(alignment.stats.output.folder, basename(stats))
  cmd <- paste(cmd, "; mv", stats, stats.new)

  # Wrap commands output ----
  cmd <- data.table(file.type= c("bam", "align.stats"),
                    path= c(bam, stats.new),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
