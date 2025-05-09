#' Import BAM File Using Rsamtools
#'
#' Imports a BAM file using Rsamtools and returns its content as a data.table.
#'
#' @param file A character string specifying the path to the BAM file.
#' @param sel A character vector specifying the columns to import. Defaults to Rsamtools::scanBamWhat().
#' @param new.col.names A character vector specifying the new column names for the imported data.
#' Defaults to `c("readID", "samFlag", "seqnames", "strand", "leftStart", "width", "mapq", "cigar",
#' "mateSeqnames", "mateLeftStart", "insertSize", "seq", "qualScore", "groupID", "mate_status")`.
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
                      new.col.names= c("readID",
                                       "samFlag",
                                       "seqnames",
                                       "strand",
                                       "leftStart",
                                       "width",
                                       "mapq",
                                       "cigar",
                                       "mateSeqnames",
                                       "mateLeftStart",
                                       "insertSize",
                                       "seq",
                                       "qualScore",
                                       "groupID",
                                       "mate_status"))
{
  # Checks ----
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
