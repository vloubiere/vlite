#' Generate Commands for Assigning Insertions to Genomic Features
#'
#' @description
#' Creates shell commands to assign insertions from a BAM file to genomic features using a GTF annotation file.
#' Outputs a BED file of unique insertions and a counts file with assigned reads.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param output.prefix Prefix for the output files. If not provided, it is derived from the input BAM filename.
#' @param genome Reference genome name (e.g., "mm10", "hg38"). If not provided, gtf must be specified.
#' @param gtf Path to the GTF annotation file. Default= NULL.
#' @param bed.output.folder Directory for the BED file of unique insertions. Default= "db/bed/".
#' @param counts.output.folder Directory for the counts file. Default= "db/counts/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file labels ("bed.file", "counts.file").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run the insertion assignment pipeline.
#'
#' @examples
#' # Assign insertions for a BAM file using the mm10 genome
#' cmd <- cmd_assignInsertions(
#'   bam = "/data/bam/sample.bam",
#'   genome = "mm10"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Assign insertions using a custom GTF file
#' cmd <- cmd_assignInsertions(
#'   bam = "/data/bam/sample.bam",
#'   gtf = "/path/to/custom.gtf"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_assignInsertions <- function(bam,
                                 output.prefix= NULL,
                                 genome,
                                 gtf,
                                 bed.output.folder= "db/bed/",
                                 counts.output.folder= "db/counts/",
                                 Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(missing(genome) & is.null(gtf))
    stop("genome is missing and gtf file was set to NULL.")

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
  counts.file.prefix <- file.path(counts.output.folder, paste0(output.prefix, "_assigned_counts"))
  fw.counts.file <- paste0(counts.file.prefix, "_same_strand.txt")
  rev.counts.file <- paste0(counts.file.prefix, "_rev_strand.txt")

  # Command ----
  cmd <- paste(Rpath,
               system.file("Rscript", "assign_ORFtag_insertions.R", package = "vlite"),
               bam,
               gtf,
               bed.file,
               counts.file.prefix)

  # Wrap commands output ----
  cmd <- data.table(file.type= c("bed.file", "fw.counts.file", "rev.counts.file"),
                    path= c(bed.file, fw.counts.file, rev.counts.file),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
