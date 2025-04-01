#' Export Genomic Ranges to BED Format Files
#'
#' @description
#' Exports genomic coordinates to standard BED, narrowPeak, or broadPeak files.
#' Automatically handles coordinate conversion and column formatting according to
#' the specified output format.
#'
#' @param bed Input genomic ranges in one of these formats:
#' \itemize{
#'   \item Character vector of genomic coordinates ("chr:start-end[:strand]")
#'   \item GRanges object
#'   \item data.frame/data.table with required columns
#'   \item Path(s) to existing .bed, .narrowPeak, or .broadPeak files
#' }
#' @param file Character. Output file path ending in .bed, .narrowPeak, or .broadPeak
#'
#' @details
#' Output format requirements:
#'
#' 1. BED format (.bed):
#'    - Required columns: seqnames, start
#'    - Auto-generated if missing:
#'      * end (set to start position)
#'      * name (set to ".")
#'      * score (set to NA)
#'      * strand (set to ".")
#'
#' 2. narrowPeak format (.narrowPeak):
#'    - Required columns:
#'      * seqnames (chr)
#'      * start (integer)
#'      * end (integer)
#'      * name (character)
#'      * score (numeric)
#'      * strand (character)
#'      * signalValue (numeric)
#'      * pValue (numeric)
#'      * qValue (numeric)
#'      * peak (integer)
#'
#' 3. broadPeak format (.broadPeak):
#'    - Same as narrowPeak but without the 'peak' column
#'
#' Note: Start coordinates are automatically converted to 0-based format
#' as required by BED specification.
#'
#' @return Writes a tab-delimited file to the specified path. No return value.
#'
#' @examples
#' # Export simple BED format
#' bed <- data.table(seqnames = "chr2L",
#'                   start = 1000,
#'                   end = 2000,
#'                   strand = "+")
#' exportBed(bed, file = "test.bed")
#' importBed("test.bed")[]
#'
#' # Export from character coordinates
#' exportBed("chr3R:1000-2000:+", file = "test.bed")
#' importBed("test.bed")[]
#'
#' # Export narrowPeak format (requires all columns)
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
  # Check output type
  type <- gsub(".*[.](.*)$", "\\1", file)

  # Import file
  current <- importBed(bed)

  # Format required columns
  if(type=="bed")
  {
    # Check format
    if(!all(c("seqnames", "start") %in% names(current)))
      stop("seqnames and start columns are required!")
    # Add missing columns
    if(!"end" %in% names(current)) current[, end := as.integer(start)]
    if(!"name" %in% names(current)) current[, name := as.character(".")]
    if(!"score" %in% names(current)) current[, score := as.numeric(NA)]
    if(!"strand" %in% names(current)) current[, strand := as.character(".")]
    # Order columns
    col.names <- c("seqnames", "start", "end", "name", "score", "strand")
    setcolorder(current, col.names)
  }else if(type=="narrowPeak")
  {
    # Check all required columns exist
    col.names <- c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
    if(!all(col.names %in% names(current)))
      stop(paste("Missing columns for narrowPeak format:",
                 paste0(setdiff(col.names, names(current)), collapse = ", ")))
    # Check columns classes
    current[, seqnames:= as.character(seqnames)]
    current[, start:= as.integer(start)]
    current[, end:= as.integer(end)]
    current[, name:= as.character(name)]
    current[, score:= as.numeric(score)]
    current[, strand:= as.character(strand)]
    current[, signalValue:= as.numeric(signalValue)]
    current[, pValue:= as.numeric(pValue)]
    current[, qValue:= as.numeric(qValue)]
    current[, peak:= as.integer(peak)]
    # Order
    current <- current[, (col.names), with= FALSE]
  }else if(type=="broadPeak")
  {
    # Check all required columns exist
    col.names <- c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
    if(!all(col.names %in% names(current)))
      stop(paste("Missing columns for braodPeak format:",
                 paste0(setdiff(col.names, names(current)), collapse = ", ")))
    # Check columns classes
    current[, seqnames:= as.character(seqnames)]
    current[, start:= as.integer(start)]
    current[, end:= as.integer(end)]
    current[, name:= as.character(name)]
    current[, score:= as.numeric(score)]
    current[, strand:= as.character(strand)]
    current[, signalValue:= as.numeric(signalValue)]
    current[, pValue:= as.numeric(pValue)]
    current[, qValue:= as.numeric(qValue)]
    # Order
    current <- current[, (col.names), with= FALSE]
  }else
    stop("Path should have .bed, .narrowPeak or .broadPeak extension")

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
