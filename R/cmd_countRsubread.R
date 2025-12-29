#' Generate Commands for RNA-Seq Read Counting Using Rsubread
#'
#' @description
#' Creates shell commands to count RNA-Seq reads from a BAM file using the Rsubread ?Rsubread::featureCounts() function.
#' Outputs a counts file and a statistics file.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param layout Sequencing layout, either "SINGLE" or "PAIRED".
#' @param output.prefix Prefix for the output files. If not provided, it is derived from the input BAM filename.
#' @param genome Reference genome name (e.g., "mm10", "dm6"). If not provided, an annotation.file must be specified.
#' @param annotation.file Path to a '.gtf', '.gtf.gz' or '.saf' annotation file. If genome is specified,
#' the corresponding gtf file will be used. SAF files are 1-based (unlike bed) and should include the following columns: GeneID, Chr, Start, End, Strand. Default= NULL.
#' @param GTF.attrType.extra Additional GTF attribute to include in the output (e.g., gene symbol). Default= NULL.
#' @param allowMultiOverlap If a read overlaps more than one feature, should be assigned to all overlapping features? Default= FALSE.
#' @param counts.stats.output.folder Directory for the alignment statistics file. Default= "db/count_stats/".
#' @param counts.output.folder Directory for the counts file. Default= "db/counts/".
#' @param Rpath Path to the Rscript binary. Default= "Rscript".
#' @param cores Number of CPU cores to use. Default= 8.
#'
#' @return A data.table with:
#' - `file.type`: Output file labels ("count.stats", "counts").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run the RNA-Seq read counting pipeline.
#' - `job.name`: Default name for the job = "countRsubReads".
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
#'   annotation.file = "/path/to/custom.gtf"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_countRsubread <- function(bam,
                              layout,
                              output.prefix= NULL,
                              genome,
                              annotation.file= NULL,
                              GTF.attrType.extra= NULL,
                              allowMultiOverlap= FALSE,
                              counts.stats.output.folder= "db/count_stats/",
                              counts.output.folder= "db/counts/",
                              Rpath= "Rscript",
                              cores= 4)
{
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided. See ?Rsubread::featureCounts() for alternatives.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(missing(genome)) {
    if(is.null(annotation.file)) {
      stop("genome is missing and no gtf/saf annotation.file is provided.")
    } else if(!grepl(".gtf$|.gtf.gz$|.saf$", annotation.file)) {
      stop("annotation file should be in '.gtf', '.gtf.gz' or '.saf' format.")
    }
  }

  # Retrieve gtf ----
  if(!missing(genome)) {
    gtf <- switch(
      genome,
      "mm10"= "/groups/stark/vloubiere/genomes/Mus_musculus/GENCODE/gencode.vM25.protein_coding_mRNASeq.gtf.gz",
      "dm6"= "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/flybase/dm6/dmel-mRNA-r6.36.gtf"
    )
    # Extra column to be added to gene ID (typically, gene symbol)
    GTF.attrType.extra <- switch(
      genome,
      "mm10"= "gene_name",
      "dm6"= "gene_symbol"
    )
  } else {
    genome <- sub("\\.[^.]*$", "", annotation.file)
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
    annotation.file,
    stats.file,
    counts.file,
    allowMultiOverlap,
    GTF.attrType.extra
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("count.stats", "counts"),
                    path= c(stats.file, counts.file),
                    cmd= cmd,
                    cores= cores,
                    job.name= "countRsubReads")

  # Return ----
  return(cmd)
}
