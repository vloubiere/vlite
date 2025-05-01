#' Export Genomic Ranges to BED Format Files
#'
#' @description
#' Exports genomic ranges as .bed, .narrowPeak or .broadPeak files
#' File specifications can be found at https://genome.ucsc.edu/FAQ/FAQformat.html.
#'
#' @param bed Input genomic ranges in one of the following formats:
#' \itemize{
#'   \item data.table or data.frame with required columns.
#'   \item Character vector of genomic coordinates ("chr:start-end[:strand]").
#'   \item GRanges object.
#'   \item File path(s) to .bed, .narrowPeak, or .broadPeak files. If several files are provided, they will be catenated.
#' }
#' @param file Output file path ending with .bed, .narrowPeak or .broadPeak extensions.
#'
#' @details
#' Upon saving, 1-based start coordinates will be converted to 0-based, following BED format specifications.
#'
#' @examples
#' # Export simple BED format:
#' bed <- data.table(seqnames = "chr2L",
#'                   start = 1000,
#'                   end = 2000,
#'                   strand = "+")
#' exportBed(bed, file = "test.bed")
#' exportBed("chr3R:1000-2000:+", file = "test.bed")
#'
#' # Export narrowPeak format:
#' peak_data <- data.table(
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
#' exportBed(peak_data, file = "test.narrowPeak")
#'
#' @export
exportBed <- function(bed, file)
{
  # Output type
  type <- gsub(".*[.](.*)$", "\\1", file)
  if(!type %in% c("bed", "narrowPeak", "broadPeak"))
    stop("Path should have .bed, .narrowPeak or .broadPeak extension.")

  # Import file (also checks columns classes are correct)
  current <- importBed(bed)
  
  # Add missing columns for bed format
  if(type=="bed") {
    if(!"end" %in% names(current)) current[, end := as.integer(start)]
    if(!"name" %in% names(current)) current[, name := as.character(".")]
    if(!"score" %in% names(current)) current[, score := as.numeric(NA)]
    if(!"strand" %in% names(current)) current[, strand := as.character(".")]
  }

  # Check if required columns exist
  col.names <- if(type=="bed") {
    c("seqnames", "start", "end", "name", "score", "strand")
  } else if(type=="broadPeak") {
    c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
  } else if(type=="narrowPeak") {
    c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
  }
  if(!all(col.names %in% names(current))) {
    missing.cols <- setdiff(col.names, names(current))
    missing.cols <- paste0(missing.cols, collapse = ", ")
    stop(paste("Missing columns for" , type, "format:", missing.cols))
  }
  
  # Select columns and order them
  current <- current[, (col.names), with= FALSE]

  # Convert start to 0 base and make sure coor are integers
  current[, start:= start-1]

  # Set scipen to avoid scientific notation
  options(scipen = 999)

  # Save
  fwrite(current,
         file = file,
         col.names = FALSE,
         row.names = FALSE,
         sep= "\t",
         quote= FALSE,
         na= NA)
}
