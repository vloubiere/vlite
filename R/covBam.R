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
#' @param output.file.txt Path to a file with a .txt extension where output will be saved.
#' If this is specified, the job will be submitted to the server and the submitted command will
#' be returned. Useful for big files.
#' @param mem When output.file.txt is specified, memory to use for the job in G. Default= 32.
#' @param logs When output.file.txt is specified, folder where log files should be saved.
#' @param overwrite In case where output.file.txt is specified and a file with this name already
#' exists, should it be overwritten? Default= FALSE.
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
                   output.file.txt= NULL,
                   mem= 32,
                   logs= "db/logs/bamCov/",
                   overwrite= FALSE,
                   bedtools_path= "/software/2020/software/bedtools/2.27.1-foss-2018b/bin/bedtools")
{
  # Checks ----
  if(any(grepl(".bed$", bed))) {
    if(length(bed)>1)
      stop("A unique bed file path should be provided.")
  } else {
    bed <- importBed(bed)
    tmp <- tempfile(tmpdir = dirname(bam),
                    fileext = ".tmp.bed")
    rtracklayer::export(bed, tmp)
    bed <- tmp
  }
  if(!is.null(output.file) && !grepl(".txt$", output.file))
    stop("output.file should be in .txt format.")

  # Compose the bedtools command ---
  cmd <- paste(bedtools_path,
               "coverage -counts",
               ifelse(!ignore.strand, "-s", ""),
               "-a", bed, "-b", bam,
               ifelse(is.null(output.file.txt), "", paste(">", output.file.txt)))

  # If output file is specified ----
  if(!is.null(output.file.txt)) {
    # Command object
    cmd <- data.table(file.type= "counts.file",
                      path= output.file.txt,
                      cmd= cmd,
                      cores= 1,
                      job.name= "covBam")
    # Submit
    vl_submit(cmd,
              mem= mem,
              logs = logs,
              overwrite = overwrite,
              execute = TRUE)
    # Return cmd
    return(cmd)
  } else {
    # Compute and import tmp output ----
    cov <- data.table::fread(cmd = cmd, header = FALSE)[[7]] # Counts stored in column 7
    # Remove temp file
    if(exists("tmp"))
      file.remove(tmp)
    # Return cov
    return(cov)
  }
}
