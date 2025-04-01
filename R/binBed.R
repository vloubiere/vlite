#' Bin Genomic Regions into Fixed-Size or Sliding Windows
#'
#' @description
#' Creates fixed-size or sliding window bins from genomic regions. Supports both
#' equal-division binning and sliding window approaches, with options for
#' strand-awareness and centered binning.
#'
#' @param bed Input genomic ranges in any format compatible with `importBed()`:
#' \itemize{
#'   \item Character vector of ranges ("chr:start-end[:strand]")
#'   \item GRanges object
#'   \item data.frame/data.table with required columns
#'   \item Path to a BED file
#' }
#' @param nbins Integer. Number of equal-sized bins per region. If specified,
#'   `bins.width` is ignored.
#' @param bins.width Integer. Size of each bin in base pairs. Used only when
#'   `nbins` is not specified.
#' @param steps.width Integer. Distance between starts of consecutive bins.
#'   Default equals `bins.width` (non-overlapping bins). Smaller values create
#'   overlapping bins.
#' @param bins.width.min Logical. If `TRUE`, only returns complete bins matching
#'   `bins.width` exactly. Default is `FALSE`.
#' @param centered Logical. If `TRUE`, centers bins around region midpoints.
#'   Default is `FALSE`.
#' @param ignore.strand Logical. If `TRUE`, bins left-to-right regardless of
#'   strand. Default is `FALSE`.
#'
#' @details
#' **Binning Methods**:
#'
#' 1. Equal Division (`nbins`):
#'    - Splits each region into `nbins` equal parts
#'    - Useful for comparing regions of different sizes
#'    - All bins within a region have equal width
#'
#' 2. Sliding Window (`bins.width`):
#'    - Creates fixed-width bins
#'    - Controls overlap via `steps.width`:
#'      * `steps.width = bins.width`: Non-overlapping bins
#'      * `steps.width < bins.width`: Overlapping bins
#'      * `steps.width > bins.width`: Gaps between bins
#'    - Option to filter incomplete bins (`bins.width.min`)
#'
#' **Strand Handling**:
#'
#' - `ignore.strand = FALSE`:
#'   * Positive/unstranded: Left to right binning
#'   * Negative strand: Right to left binning
#' - `ignore.strand = TRUE`: Always bins left to right
#'
#' **Centered Binning**:
#'
#' When `centered = TRUE`:
#' - Extends regions symmetrically around midpoint
#' - Useful for analyzing features relative to central points
#' - Extension size depends on binning method:
#'   * Equal division: (region_width/nbins)/2
#'   * Sliding window: bins.width/2
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item line.idx: Index linking bins to source regions
#'   \item bin.idx: Sequential bin number within each region
#'   \item seqnames: Chromosome/sequence name
#'   \item start: Bin start position
#'   \item end: Bin end position
#'   \item Additional columns from input preserved
#' }
#'
#' @examples
#' # Equal division binning
#' regions <- c(
#'   "chr2L:1-10:+",    # 10bp region
#'   "chr2L:100-200:-"  # 100bp region
#' )
#'
#' # Split into 2 equal bins
#' binBed(regions, nbins = 2)
#'
#' # Split into 5 equal bins
#' binBed(regions, nbins = 5)
#'
#' # Sliding window with different overlaps
#' region <- "chr2L:1-20:+"
#'
#' # Non-overlapping 5bp bins
#' binBed(region, bins.width = 5, steps.width = 5)
#'
#' # Overlapping bins (2bp step)
#' binBed(region, bins.width = 5, steps.width = 2)
#'
#' # With gaps (7bp step)
#' binBed(region, bins.width = 5, steps.width = 7)
#'
#' # Centered binning
#' tss <- "chr2L:1000:+"  # TSS position
#' # Create 100bp bins centered on TSS
#' binBed(tss, bins.width = 100, centered = TRUE)
#'
#' # Strand-specific binning
#' binBed("chr2L:1-10:+", nbins = 2)  # Left to right
#' binBed("chr2L:1-10:-", nbins = 2)  # Right to left
#'
#' @export
binBed <- function(bed,
                   nbins,
                   bins.width = NULL,
                   steps.width = bins.width,
                   bins.width.min= FALSE,
                   centered= FALSE,
                   ignore.strand= FALSE)
{
  # Checks
  method <- if(!missing(nbins)) {
    "nbins"
  } else if(is.numeric(bins.width) & is.numeric(steps.width)) {
    "slidingWindow"
  } else
    stop("nbins or bins.width should be specified!")

  # Hard copy for incapsulation ----
  regions <- importBed(bed)

  # Add index and negative strand columns ----
  regions[, line.idx:= .I]
  regions[, neg.strand:= FALSE]
  if(!ignore.strand && "strand" %in% names(regions))
    regions[, neg.strand:= strand=="-"]

  # Resize if bins have to be centered ----
  if(centered) {
    # Compute ext size
    ext <- if(method=="nbins") {
      floor(regions[, (end-start+1)]/nbins/2)
    } else if(method=="slidingWindow") {
      floor(bins.width/2)
    }
    # Extend
    regions <- resizeBed(regions,
                         center = "region",
                         upstream = ext,
                         downstream = ext,
                         ignore.strand = ignore.strand)
  }

  # Compute bins ----
  bins <- if(method=="nbins") {

    # Using fixed number of bins ----
    .bin_regions_using_Nbins(regions= regions,
                             nbins= nbins)

  } else if(method=="slidingWindow"){

    # Using fixed bin widths and steps ----
    .bin_regions_using_width(regions= regions,
                             bins.width= bins.width,
                             steps.width= steps.width,
                             bins.width.min= bins.width.min)
  }

  # Add line.idx and bin.idx to original bed ----
  input <- importBed(bed)
  input <- input[bins$line.idx, !c("start", "end")]
  res <- cbind(bins, input)

  # Order columns ----
  setcolorder(res,
              c("line.idx", "bin.idx", "seqnames"))

  # Return bins and a ----
  return(res)
}
