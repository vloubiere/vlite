#' Generate Commands for Counting PRO-seq reads
#'
#' @description
#' This function generates shell commands to count PRO-seq reads from a UMI count file using a .rds genomic annotation file.
#' The output includes a counts table with assigned reads for the specified genomic features.
#'
#' @param umi.count.file Path to the input UMI count file. Must be a single .txt file.
#' @param annotation.file Path to the annotation file in .rds format containing the following columns:
#' seqnames, start, end, strand, cluster.id. The mm10 annotations are stored in '/groups/stark/vloubiere/genomes/Mus_musculus/PROseq/'.
#' @param blacklist.file An optional file of blacklisted regions for which overlapping reads should be removed.
#' Default= NULL. For mm10, a list of blacklisted regions (including tRNAs/snRNAs...) is availaible in
#' '/groups/stark/vloubiere/genomes/Mus_musculus/PROseq/'.
#' @param feature Name of the genomic features in the annotation.file (promoter, body...). If not provided, it is derived from the annotation.file name.
#' @param output.prefix Prefix for the output files. If not provided, it is derived from the umi.count.file name.
#' @param count.tables.output.folder Directory for the output counts table. Default: "db/counts/PROseq/".
#' @param Rpath Path to the Rscript binary. Default: "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with the following columns:
#' - `file.type`: Label for the output file ("count.table").
#' - `path`: Path to the output counts table.
#' - `cmd`: Shell command to run the PRO-seq read counting pipeline.
#' - `job.name`: Default name for the job = "countProseqReads".
#'
#' @examples
#' # Example usage
#' cmd <- cmd_countPROseqReads(
#'   umi.count.file = "example_umi_counts.txt",
#'   annotation.file = "/groups/stark/vloubiere/genomes/Mus_musculus/PROseq/mm10_transcript.rds",
#'   feature = "promoters",
#'   output.prefix = "example_output"
#' )
#' print(cmd)
#'
#' @export
cmd_countPROseqReads <- function(umi.count.file,
                                 annotation.file,
                                 blacklist.file= NULL,
                                 feature= NULL,
                                 output.prefix= NULL,
                                 count.tables.output.folder= "db/count_tables/PROseq/",
                                 Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if umi.count.file exists to allow wrapping!) ----
  if(length(umi.count.file)!=1)
    stop("A unique umi.count.file should be provided.")
  if(!grepl(".txt$", umi.count.file))
    stop("umi.count.file should be in .txt format")
  if(is.null(output.prefix))
    output.prefix <- gsub(".txt$", "", basename(umi.count.file))
  if(length(annotation.file)>1)
    stop("A single annotation file should be provided.")
  if(!grepl(".rds$", annotation.file))
    stop("The annotation file should be in .rds format.")
  check.annot <- readRDS(annotation.file)
  if(!is.data.table(check.annot))
    stop("annotation.file should contain a data.table")
  if(!all(c("seqnames", "start", "end", "strand", "cluster.id") %in% names(check.annot)))
    stop("Annotation file should contain columns 'seqnames', 'start', 'end', 'strand', 'cluster.id'")
  if(!is.null(blacklist.file)) {
    if(!grepl(".rds$", blacklist.file))
      stop("The blacklist file should be in .rds format.")
    check.blacklist <- readRDS(blacklist.file)
    if(!is.data.table(check.blacklist))
      stop("blacklist.file should contain a data.table")
    if(!all(c("seqnames", "start", "end", "strand", "cluster.id") %in% names(check.blacklist)))
      stop("Blacklist file should contain columns 'seqnames', 'start', 'end', 'strand', 'cluster.id'")
  }
  if(is.null(feature))
    feature <- gsub(".rds$", "", basename(annotation.file))

  # Output files ----
  count.table <- file.path(count.tables.output.folder, paste0(output.prefix, "_", feature, "_counts.txt"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "count_PROseq_reads.R", package = "vlite"),
    umi.count.file,
    annotation.file,
    count.table,
    blacklist.file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "count.table",
                    path= count.table,
                    cmd= cmd,
                    job.name= "countProseqReads")
  if(!is.null(blacklist.file)) {
    blacklist.counts.file <- gsub(".txt$", "__blacklisted.txt", count.table)
    cmd <- rbind(cmd, cmd)
    cmd[2, c("file.type", "path"):= .("blacklist.counts", blacklist.counts.file)]
  }

  # Return ----
  return(cmd)
}
