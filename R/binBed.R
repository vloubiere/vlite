#' Bin Genomic Regions into Fixed-Size or Sliding Windows
#'
#' @description
#' Creates fixed-size or sliding window bins from genomic ranges. Supports both
#' equal-division binning and sliding window approaches.
#'
#' @param bed Input genomic ranges in any format compatible with ?importBed().
#' @param nbins Integer specifying the number of equal-sized bins in which each region is divided.
#' @param center.nbins If set to TRUE, bins are centered around the regions' midpoints. If set to
#' FALSE (default), the first bin will start from the most upstream coordinate.
#' Only meaningful when nbins is specified.
#' @param bins.width Integer specifying the width of the bins in which each region is divided. 
#' Used only when nbins is not specified. See bins.width.min for handling shorter bins
#' near bed boundaries. 
#' @param steps.width Integer specifying the distance between starts of consecutive bins. Used in
#' combination with bins.width. Default= bins.width, resulting in non-overlapping bins.
#' Smaller values create overlapping bins.
#' If specified, bins.width is ignored.
#' @param bins.width.min If set to TRUE, only bins matching the size specified in bins.width
#' are returned. This is useful to exclude shorter bins nearby regions' boundaries.
#' Default= FALSE.
#' @param ignore.strand If set to TRUE, bins are oriented in 5' to 3' positions, regardless of 
#' their strand. If set to FALSE (default), bins orientation follow input regions' strand.
#'
#' @details
#'
#' **Centered Binning**:
#' Centering bins (center.nbins = TRUE) is useful for analyzing features relative to central points.
#' In this case, each region will first be extended symmetrically by region_width/(nbins-1)/2.
#' Only available when using nbins.
#'
#' @return A gr data.table with columns:
#' \itemize{
#'   \item line.idx: line index of the region in bed.
#'   \item bin.idx: bin index (unique for each line.idx).
#'   \item seqnames: chromosome or sequence name.
#'   \item start: bin start position.
#'   \item end: bin end position.
#'   \item Additional columns within bed input are preserved.
#' }
#'
#' @examples
#' # Equal division binning
#' regions <- importBed(c("chr2L:1-10:+", "chr2L:100-200:-"))
#'
#' # Split into 5 equal bins
#' binBed(regions, nbins = 5)
#'
#' # Non-overlapping 50bp bins
#' binBed(regions, bins.width = 50, steps.width = 50)
#'
#' # Overlapping 50bp bins
#' binBed(regions, bins.width = 50, steps.width = 25)
#'
#' # Compare start-anchored vs. centered binning
#' nc1 <- binBed(regions, nbins= 5)
#' c1 <- binBed(regions, nbins= 5, center.nbins = TRUE)
#' nc2 <- binBed(regions, bins.width = 25)
#' 
#' plot(x= c(50, 250), y= c(4, 10), type= "n", xlab= "coordinates", ylab= NA, yaxt= "n")
#' abline(v= c(100, 200))
#' nc1[line.idx==2][, {text(x= mean(c(start[1], end[.N])), y= 9, pos= 3, "nbins=5"); rect(start, 8, end, 9)}]
#' c1[line.idx==2][, {text(x= mean(c(start[1], end[.N])), y= 7, pos= 3, "centered"); rect(start, 6, end, 7)}]
#' nc2[line.idx==2][, {text(x= mean(c(start[1], end[.N])), y= 5, pos= 3, "bins.width=25"); rect(start, 4, end, 5)}]
#'
#' @export
binBed <- function(bed,
                   nbins,
                   center.nbins= FALSE,
                   bins.width = NULL,
                   steps.width = bins.width,
                   bins.width.min= FALSE,
                   genome= NULL,
                   ignore.strand= FALSE)
{
  # Checks ----
  method <- if(!missing(nbins)) {
    if(nbins %% 1 != 0 | nbins<1) {
      stop("nbins should be an integer >= 1!")
    } else
      "nbins"
  } else if(is.numeric(bins.width) & is.numeric(steps.width)) {
    if(any(c(bins.width, steps.width) %% 1 != 0) || any(c(bins.width, steps.width)<1)) {
      stop("bins.width and steps.width should be integers >= 1!")
    } else
      "slidingWindow"
  } else
    stop("nbins or bins.width should be specified!")
  if(center.nbins && method!="nbins")
    stop("The 'center.nbins' option is only meaningful when 'nbins' is provided.")

  # Hard copy for incapsulation ----
  regions <- importBed(bed)

  # Add index and negative strand columns ----
  regions[, line.idx:= .I]
  regions[, neg.strand:= FALSE]
  if(!ignore.strand && "strand" %in% names(regions))
    regions[, neg.strand:= strand=="-"]

  # Resize if bins have to be centered ----
  if(center.nbins) {
    # Compute ext size
    ext <- floor(regions[, (end-start+1)]/(nbins-1)/2)
    # Extend
    regions <- resizeBed(regions,
                         center = "region",
                         upstream = ext,
                         downstream = ext,
                         genome= genome,
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
