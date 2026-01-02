#' Import Genomic Ranges in BED Format
#'
#' @description
#' Imports genomic coordinates from various input formats into a standardized data.table format.
#' Supported input types: GRanges objects, character strings ("chr2L:30000-31000:-"), .bed, .narrowPeak
#' or .broadPeak file paths (see details).
#'
#' @param bed Input genomic ranges in one of the following formats:
#' \itemize{
#'   \item Character vector of genomic coordinates ("chr:start-end[:strand]").
#'   \item File path to .bed, .narrowPeak, or .broadPeak file.
#'   \item GRanges object.
#'   \item data.frame or data.table with required columns.
#' }
#'
#' @details
#' Character coordinates should follow the syntax "chr:start-end[:strand]".
#' File specifications can be found at https://genome.ucsc.edu/FAQ/FAQformat.html.
#' GRanges objects are directly converted to data.table and metadata columns are preserved.
#' Required columns for data.frames/data.tables are 'seqnames', 'start', 'end' and 'strand'.
#' If missing, end and strand will be set to 'start' and '*', respectively.
#'
#' @examples
#' # Import from character coordinates
#' importBed(c("chr2L:1000-2000:+", "chr3L:1000-2000:-", "chr3L:4000-6000:*"))[]
#'
#' # From GRanges
#' gr <- importBed(GRanges("chr2L", IRanges(c(1000, 5000), c(2000, 6000))))[]
#'
#' # From bed file
#' bed_file <- tempfile(fileext = ".bed")
#' exportBed(gr, bed_file)
#' importBed(bed_file)[]
#'
#' # From narrowpeak file
#' peaks <- data.table(
#'   seqnames = "chr2L",
#'   start = 1000,
#'   end = 2000,
#'   name = "peak1",
#'   score = 100,
#'   strand = "+",
#'   signalValue = 5.5,
#'   pValue = 0.001,
#'   qValue = 0.05,
#'   peak = 1500
#' )
#' narrowPeak_file <- tempfile(fileext = ".narrowPeak")
#' exportBed(peaks, narrowPeak_file)
#' importBed(narrowPeak_file)
#'
#' @export
importBed <- function(bed)
{
  if(is.character(bed)) {
    if(any(grepl("\\.(bed|narrowPeak|broadPeak|bed.gz|narrowPeak.gz|broadPeak.gz)$", bed))) {

      # Import bed file using rtracklayer ----
      current <- data.table::as.data.table(rtracklayer::import(bed))

    } else {

      # Import character coordinates using GenomicRanges ----
      gr <- gsub(",", "", bed) # Remove potential comas
      gr <- GenomicRanges::GRanges(gr)
      check.col <- names(GenomicRanges::mcols(gr))
      current <- data.table::as.data.table(gr)
      # If no 'width' column in mcols, remove the one created during data.table conversion
      if(!length(check.col) || !"width" %in% check.col)
        current$width <- NULL
      # Store the original coordinates as names
      current$name <- bed
    }
  } else if(is.data.frame(bed)) {

    # Coerce to data.table ----
    current <- if(data.table::is.data.table(bed)) data.table::copy(bed) else data.table::as.data.table(bed)
    # Check required columns
    if(!all(c("seqnames", "start") %in% names(current)))
      stop("Genomic ranges should contain at least 'seqnames' and 'start' fields.")
    # Format required columns
    current[, seqnames := as.character(seqnames)]
    current[, start := as.integer(start)]
    if(!"end" %in% names(current)) current[, end:= start]
    current[, end:= as.integer(end)]
    if(!"strand" %in% names(current)) current[, strand:= "*"]
    # Check strand values
    if(any(!current$strand %in% c("+", "-", "*")))
      stop("All strand values should be one of c('+', '-' or '*') -> malformed genomic ranges!")

  } else if(inherits(bed, "GRanges")) {

    # Coerce to data.table ----
    gr <- GenomicRanges::GRanges(bed)
    check.col <- names(GenomicRanges::mcols(gr))
    current <- data.table::as.data.table(gr)
    # If no 'width' column in mcols, remove the one created during data.table conversion
    if(!length(check.col) || !"width" %in% check.col)
      current$width <- NULL

  } else {
    stop("Input format could not be determined. See ?vlite::importBed.")
  }

  # Sanity check ----
  if(any(current$start<1 | current[,start>end]))
    warning("Some regions with start<1 or start>end!")
  
  # Return bed data.table ----
  return(current)
}
