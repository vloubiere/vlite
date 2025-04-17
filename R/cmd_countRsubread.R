#' Generate Commands for RNA-Seq Read Counting Using Rsubread
#'
#' @description
#' Creates shell commands to count RNA-Seq reads from a BAM file using the Rsubread ?Rsubread::featureCounts() function.
#' Outputs a counts file and a statistics file.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param layout Sequencing layout, either "SINGLE" or "PAIRED".
#' @param output.prefix Prefix for the output files. If not provided, it is derived from the input BAM filename.
#' @param genome Reference genome name (e.g., "mm10", "dm6"). If not provided, gtf must be specified.
#' @param gtf Path to the GTF annotation file. Default= NULL.
#' @param GTF.attrType.extra Additional GTF attribute to include in the output (e.g., gene symbol). Default= NULL.
#' @param counts.stats.output.folder Directory for the alignment statistics file. Default= "db/count_stats/".
#' @param counts.output.folder Directory for the counts file. Default= "db/counts/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file labels ("count.stats", "counts").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run the RNA-Seq read counting pipeline.
#'
#' @examples
#' # Count reads for a single-end RNA-Seq BAM file
#' cmd <- cmd_countRsubread(
#'   bam = "/data/bam/sample.bam",
#'   layout = "SINGLE",
#'   genome = "mm10"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Count reads using a custom GTF file
#' cmd <- cmd_countRsubread(
#'   bam = "/data/bam/sample.bam",
#'   layout = "PAIRED",
#'   gtf = "/path/to/custom.gtf"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
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
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided. See ?Rsubread::featureCounts() for alternatives.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(missing(genome) & is.null(gtf))
    stop("genome is missing and gtf file was set to NULL.")

  # Retrieve gtf ----
  if(!missing(genome)) {
    gtf <- switch(genome,
                  "mm10"= "/groups/stark/vloubiere/genomes/Mus_musculus/GENCODE/gencode.vM25.protein_coding_mRNASeq.gtf.gz",
                  "dm6"= "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/flybase/dm6/dmel-mRNA-r6.36.gtf")
    # Extra column to be added to gene ID (typically, gene symbol)
    GTF.attrType.extra <- switch(genome,
                                 "mm10"= "gene_name",
                                 "dm6"= "gene_symbol")
  } else {
    genome <- gsub(".gtf$", "", basename(gtf))
  }

  # Output files paths ----
  stats.file <- file.path(counts.stats.output.folder, paste0(output.prefix, "_statistics.txt"))
  counts.file <- file.path(counts.output.folder, paste0(output.prefix, "_counts.txt"))

  # Count command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "count_Rsubread.R", package = "vlite"),
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
