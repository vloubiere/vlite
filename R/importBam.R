#' Import BAM File Using Rsamtools
#'
#' Imports a BAM file using Rsamtools and returns its content as a data.table.
#'
#' @param file A character string specifying the path to the BAM file.
#' @param sel A character vector specifying the columns to import. Defaults to Rsamtools::scanBamWhat().
#' @param new.col.names A character vector specifying the new column names for the imported data.
#' If set to NULL (default), 'sel' will be matched to: "readID", "samFlag", "seqnames", "strand", "leftStart",
#' "width", "mapq", "cigar", "mateSeqnames", "mateLeftStart", "insertSize", "seq", "qualScore", "groupID",
#' "mate_status".
#'
#' @return
#' A data.table containing the content of the BAM file with the specified columns and column names.
#'
#' @examples
#' # Import a BAM file with all default columns
#' importBam(file= "path/to_bam_file")
#'
#' # Import specific columns columns
#' importBam(file= "path/to_bam_file", sel = c("qname", "flag"), new.col.names = c("readID", "samFlag"))
#'
#' @export
importBam <- function(file,
                      sel= Rsamtools::scanBamWhat(),
                      new.col.names= NULL)
{
  # Checks
  if(length(file)!=1)
    stop("file should be unique.")

  # Default col.names ----
  if(is.null(new.col.names)) {
    new.col.names <- sapply(sel, function(x) {
      switch(x,
             "qname"= "readID",
             "flag"= "samFlag",
             "rname"= "seqnames",
             "strand"= "strand",
             "pos"= "leftStart",
             "qwidth"= "width",
             "mapq"= "mapq",
             "cigar"= "cigar",
             "mrnm"= "mateSeqnames",
             "mpos"= "mateLeftStart",
             "isize"= "insertSize",
             "seq"= "seq",
             "qual"= "qualScore",
             "groupid"= "groupID",
             "mate_status"= "mate_status")
    })
  }
  if(length(sel) != length(new.col.names))
    stop("'sel' and 'new.col.names' should have the same length.")

  # Import file ----
  param <- Rsamtools::ScanBamParam(what= sel)
  .c <- Rsamtools::scanBam(file, param = param)[[1]]

  # Coerce special elements (typically DNAStringSet / PhredQuality) to character ----
  toCoerce <- which(!sapply(.c, function(x) is.vector(x) | is.factor(x)))
  for(i in toCoerce){
    .c[[i]] <- as.character(.c[[i]])
  }

  # To data.table ----
  .c <- as.data.table(.c)
  setcolorder(.c, intersect(sel, names(.c)))
  setnames(.c,
           old = sel,
           new= new.col.names,
           skip_absent = TRUE)

  # Return ----
  return(.c)
}

#' Import Raw BAM File Using Samtools and fread
#'
#' Uses samtools to import a BAM file and reads the output into R using fread.
#'
#' @param file A character string specifying the path to the BAM file.
#' @param extra_arg Optional extra arguments that will be passed to `samtools view`.
#' @param headN An integer specifying the number of starting lines to import. By default, all reads will be processed.
#' @param samtools_exe_path A character string specifying the path to the samtools executable.
#' Default= "/software/2020/software/samtools/1.9-foss-2018b/bin/samtools".
#'
#' @return
#' A data table containing the imported BAM file.
#'
#' @examples
#' # Import a BAM file
#' importBamRaw("path/to/bam/file.bam")
#'
#' # Import with additional arguments
#' importBamRaw("path/to/bam/file.bam", extra_arg = "-f 4", headN = 100)
#'
#' @export
importBamRaw <- function(file,
                         extra_arg,
                         headN,
                         samtools_exe_path= "/software/2020/software/samtools/1.9-foss-2018b/bin/samtools")
{
  # Initiate command
  cmd <- paste(samtools_exe_path, "view")

  # Add extra arguments
  if(!missing(extra_arg))
    cmd <- paste(cmd,
                 extra_arg)
  cmd <- paste(cmd, file)

  # Add head if specified
  if(!missing(headN))
  {
    if(!is.integer(headN))
      headN <- as.integer(headN)
    cmd <- paste(cmd, "| head -n", headN)
  }

  # Import
  fread(cmd= cmd,
        fill= T,
        header= F,
        sep= "\t")
}
