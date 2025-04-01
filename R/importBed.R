#' Import Genomic Ranges in BED Format
#'
#' @description
#' A versatile function that imports genomic coordinates from various input formats
#' into a standardized data.table. Supports multiple input types including BED files,
#' narrowPeak files, GRanges objects, and string coordinates.
#'
#' @param bed Input genomic ranges in one of these formats:
#' \itemize{
#'   \item Character vector of genomic coordinates ("chr:start-end[:strand]")
#'   \item File path(s) to .bed, .narrowPeak, or .broadPeak files
#'   \item GRanges object
#'   \item data.frame or data.table with required columns
#' }
#' @param extra.columns Character vector specifying names for additional columns
#'   in BED files beyond the standard fields. Default is NULL.
#'
#' @details
#' Input format handling:
#'
#' 1. BED files (.bed):
#'    - Must contain at least chromosome, start, and end columns
#'    - Additional columns can be named using extra.columns
#'
#' 2. Peak files (.narrowPeak, .broadPeak):
#'    - Automatically recognizes ENCODE peak file formats
#'    - Preserves all standard peak columns
#'
#' 3. Character coordinates:
#'    - Format: "chr:start-end[:strand]"
#'    - Strand is optional (defaults to "*" if not specified)
#'    - Example: "chr2L:1000-2000:+" or "chr2L:1000-2000"
#'
#' 4. GRanges objects:
#'    - Converts directly to data.table
#'    - Preserves metadata columns
#'
#' 5. data.frame/data.table:
#'    - Must contain columns: seqnames, start, end
#'    - Strand column is optional
#'
#' @return A data.table with standardized columns:
#' \itemize{
#'   \item seqnames: Chromosome or sequence name (character)
#'   \item start: Start position (integer)
#'   \item end: End position (integer)
#'   \item strand: Strand information (character: "+", "-", or "*")
#'   \item Additional columns from input (if any)
#' }
#'
#' @examples
#' # Import from character coordinates
#' importBed(c("chr2L:1000-2000:+", "chr3L:1000-2000:-", "chr3L:4000-6000:*"))
#'
#' # Import from GRanges
#' gr <- GRanges("chr2L", IRanges(c(1000, 5000), c(2000, 6000)))
#' importBed(gr)
#'
#' # Import single BED file
#' bed_file <- system.file("extdata", "Drosophila_transcripts_r6.36.bed",
#'                        package = "dtBedTools")
#' bed <- importBed(bed_file)
#'
#' # Import and merge multiple peak files
#' peaks1 <- system.file("extdata", "ATAC_control_eye_disc.narrowPeak",
#'                      package = "dtBedTools")
#' peaks2 <- system.file("extdata", "ATAC_tumor.narrowPeak",
#'                      package = "dtBedTools")
#' combined_peaks <- importBed(c(peaks1, peaks2))
#'
#' @export
importBed <- function(bed, extra.columns= NULL)
{
  # Import as data table ----
  current <- if(is.character(bed))
  {
    if(all(grepl(".bed$", bed)))
    {
      .importBedFile(file = bed, extra.columns = extra.columns)
    }else if(all(grepl(".narrowPeak", bed)))
    {
      .importNarrowPeakFile(file = bed)
    }else if(all(grepl(".broadPeak", bed)))
    {
      .importBroadPeakFile(file = bed)
    }else if(all(grepl("^.*:.*$", bed)))
    {
      .importCharacterCoor(coor= bed)
    }else
    {
      stop("Format could not be determined.
           If you are using character coordinates,
           use the following syntax: 'chrX:1000-2000' or 'chrX:1000-2000'")
    }
  }else if(class(bed)[1]=="GRanges")
  {
    as.data.table(bed)
  }else if(is.data.table(bed))
  {
    data.table::copy(bed)
  }else if(is.data.frame(bed))
  {
    as.data.table(bed)
  }else
    stop("Input bed(s) are inconsistent. See ?importBed().")

  # Make sure that columns have expected classes ----
  current <- .checkDataTableBedColClasses(current)

  # Return bed data.table ----
  return(current)
}
