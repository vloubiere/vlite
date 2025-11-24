#' Export Genomic Ranges to BED Format Files
#'
#' @description
#' Exports genomic ranges as .bed, .narrowPeak or .broadPeak files. See details.
#'
#' @param bed Input genomic ranges, in any format compatible with ?importBed.
#' @param file Output file path ending with .bed, .narrowPeak or .broadPeak extension.
#'
#' @details
#' File specifications can be found at https://genome.ucsc.edu/FAQ/FAQformat.html.
#' Upon saving, 1-based start coordinates will be converted to 0-based, following BED format specifications.
#'
#' @examples
#' # Export simple BED format:
#' exportBed("chr2L:1000-2000:+", file = "test.bed")
#' # Lookup saved file
#' fread("test.bed")
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
#' exportBed(bed= peak_data, file = "test.narrowPeak")
#' # Lookup saved file
#' fread("test.narrowPeak")
#'
#' @export
exportBed <- function(bed, file)
{
  # Output type ----
  type <- gsub(".*[.](.*)$", "\\1", file)
  if(!type %in% c("bed", "narrowPeak", "broadPeak"))
    stop("Path should have .bed, .narrowPeak or .broadPeak extension.")

  if(type=="bed") {

    # Export using rtracklayer ----
    rtracklayer::export(GenomicRanges::GRanges(bed), file)

  } else {

    # Custom method for broadPeak/narrowPeak ----
    current <- vlite::importBed(bed)
    # Check if required columns exist
    col.names <- if(type=="broadPeak") {
      c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
    } else if(type=="narrowPeak") {
      c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
    }
    missing.cols <- setdiff(col.names, names(current))
    if(length(missing.cols))
      stop(paste("Missing columns for" , type, "format:", paste0(missing.cols, collapse = ", ")))
    # Select and order columns
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
}
