#' Title
#'
#' A wrapper around 'bedtools coverage' that computes, for each genomic range in a
#' bed file, the number of overlapping reads in a bam file.
#'
#' @param bed Regions in any format compatible with ?importBed.
#' @param bam Path to a bam file from which overlapping reads will be counted.
#' @param ignore.strand If set to FALSE, only the reads features that are on the same strand
#' as the bed regions will be counted. If set to TRUE (default), reads are counted irrespective
#' of their strand.
#' @param output.prefix Prefix for the output file name. If specified, a 'cmd' object compatible with
#' ?vl_submit will be returned. If set to NULL (default), the command will be submitted to the system
#' using fread(cmd= cmd).
#' @param output.folder Directory for the output count files. Default: "db/counts/".
#' @param bedtools_path Path to the bedtools bin folder.
#' Default= "/software/2020/software/bedtools/2.27.1-foss-2018b/bin/bedtools".
#'
#' @return A vector of integer counts.
#' @export
#'
#' @examples
#' covBam(bed= "path/to/bed/file.bed",
#' bam = "path/to/bam/file.bam")
#'
covBam <- function(bed,
                   bam,
                   ignore.strand= TRUE,
                   output.prefix= NULL,
                   output.folder= "db/counts/",
                   bedtools_path= "/software/2020/software/bedtools/2.27.1-foss-2018b/bin/bedtools")
{
  # Checks ----
  if(any(grepl(".bed$", bed))) {
    if(length(bed)>1)
      stop("A unique bed file path should be provided.")
  } else {
    # If bed is is not a file, save it as tmp file ----
    bed <- importBed(bed)
    tmp <- tempfile(tmpdir = dirname(bam),
                    fileext = ".tmp.bed")
    rtracklayer::export(bed, tmp)
    bed <- tmp
  }

  # Output file path ----
  output.file <- if(!is.null(output.prefix))
    file.path(output.folder, paste0(output.prefix, "_counts.txt")) else
      NULL

  # Compose the bedtools command ---
  cmd <- paste(bedtools_path,
               "coverage -counts",
               ifelse(!ignore.strand, "-s", ""),
               "-a", bed, "-b", bam,
               ifelse(is.null(output.file), "", paste(">", output.file)))

  # If output file is specified ----
  if(!is.null(output.prefix)) {
    # Create command
    cmd <- data.table(file.type= "counts.file",
                      path= output.file,
                      cmd= cmd,
                      cores= 1,
                      job.name= "covBam")
    # Return cmd
    return(cmd)
  } else {

    # Compute and import counts ----
    cov <- data.table::fread(cmd = cmd,
                             header = FALSE)[[7]] # Counts stored in column 7
    # Remove temp file
    if(exists("tmp"))
      file.remove(tmp)

    # Return cov
    return(cov)
  }
}
