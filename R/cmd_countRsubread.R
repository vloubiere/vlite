#' Title
#'
#' @param bam
#' @param layout
#' @param genome
#' @param gtf
#' @param GTF.attrType.extra
#' @param counts.stats.output.folder
#' @param counts.output.folder
#' @param Rpath
#'
#' @return
#' @export
#'
#' @examples
cmd_countRsubread <- function(bam,
                              layout,
                              output.prefix= NULL,
                              genome,
                              gtf= NULL,
                              GTF.attrType.extra= NULL,
                              counts.stats.output.folder= "db/count_stats/",
                              counts.output.folder= "db/counts/",
                              Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided. See ?Rsubread::featureCounts() for alternatives.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(missing(genome) & is.null(gtf))
    stop("genome is missing and gtf file was set to NULL.")
  if(!dir.exists(counts.stats.output.folder))
    dir.create(counts.stats.output.folder, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(counts.output.folder))
    dir.create(counts.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Retrieve gtf ----
  if(!missing(genome)) {
    gtf <- switch(genome,
                  "mm10"= "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/gencode.vM25.basic.annotation.gtf.gz",
                  "dm6"= "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/flybase/dm6/dmel-all-r6.36.gtf")
    # Extra column to be added to gene ID (typically, gene symbol)
    GTF.attrType.extra <- switch(genome,
                                 "mm10"= "gene_name",
                                 "dm6"= "gene_symbol")
  } else {
    genome <- gsub(".gtf$", "", basename(gtf))
  }

  # Output files paths ----
  stats.file <- file.path(counts.stats.output.folder, paste0(output.prefix, "_", genome, "_statistics.txt"))
  counts.file <- file.path(counts.output.folder, paste0(output.prefix, "_", genome, "_counts.txt"))

  # Count command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "count_Rsubread.R", package = "genomicsPipelines"),
    bam,
    layout,
    gtf,
    stats.file,
    counts.file,
    GTF.attrType.extra
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("count.stats", "counts"),
                    path= c(stats.file, counts.file),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
