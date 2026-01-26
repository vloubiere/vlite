#' Bin Genomic Regions into Fixed-Size or Sliding Windows
#'
#' @description
#' A wrapper around `GenomicRanges::tile()` and `GenomicRanges::slidingWindows()`, that creates fixed-size
#' or sliding window bins from genomic ranges.
#'
#' @param bed Input genomic ranges in any format compatible with ?importBed.
#' @param nbins Integer specifying the number of equal-sized bins in which each region is divided.
#' @param center.nbins If set to TRUE and nbins is specified, bins are centered around the regions' midpoints.
#' If set to FALSE (default), the first bin will start from the leftmost coordinate.
#' @param genome When center.bins is set to TRUE, a BS genome name (e.g. 'dm6', 'mm10') that will be used to clip
#' resized regions that might extend beyond chromosome size (see details).
#' @param bins.width Integer specifying the width of the bins in which each region is divided.
#' Used only when nbins is not specified.
#' @param steps.width Integer specifying the distance between starts of consecutive bins. Used in
#' combination with bins.width. Default= bins.width, resulting in non-overlapping bins.
#' Smaller values create overlapping bins.
#' @param ignore.strand Bins are always generated left-to-right. If ignore.strand=FALSE (default), bin.idx is
#' reversed on `-` strand so indices increase 5′→3′.
#'
#' @details
#' **Centered Binning**:
#' Centering bins (center.nbins = TRUE) is useful for analyzing features relative to central points.
#' In this case, each region will first be extended symmetrically by region_width/(nbins-1)/2, and a
#' genome can be specified to avoid regions outside of chromosome sizes.
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
#' nbins <- 5
#' bins.width <- 25
#' nc1 <- binBed(regions, nbins= nbins)
#' c1 <- binBed(regions, nbins= nbins, center.nbins = TRUE)
#' nc2 <- binBed(regions, bins.width = bins.width)
#' 
#' plot(x= c(50, 250), y= c(4, 10), type= "n", xlab= "coordinates", ylab= NA, yaxt= "n")
#' abline(v= c(100, 200))
#' nc1[line.idx==2][, {text(x= mean(c(start[1], end[.N])), y= 9, pos= 3, paste0("nbins=", nbins)); rect(start, 8, end, 9)}]
#' c1[line.idx==2][, {text(x= mean(c(start[1], end[.N])), y= 7, pos= 3, "centered"); rect(start, 6, end, 7)}]
#' nc2[line.idx==2][, {text(x= mean(c(start[1], end[.N])), y= 5, pos= 3, paste0("bins.width=", bins.width)); rect(start, 4, end, 5)}]
#'
#' @export
binBed <- function(bed,
                   nbins,
                   center.nbins= FALSE,
                   genome= NULL,
                   bins.width = NULL,
                   steps.width = bins.width,
                   ignore.strand= FALSE)
{
  # Check method ----
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
    stop("The 'center.nbins' option only available when 'nbins' is provided.")
  
  # Import bed for incapsulation ----
  bed <- vlite::importBed(bed)
  if(any(c("line.idx", "bin.idx") %in% names(bed)))
    warning("'line.idx' or 'bin.idx' column(s) already exist in the input bed and will be overwritten.")
  
  # Add line index ----
  bed[, line.idx:= .I]
  
  # Resize if bins have to be centered ----
  if(center.nbins && nbins>1) {
    # Compute ext size
    ext <- bed[, ((end-start+1)/(nbins-1)) %/% 2L]
    # Extend
    bed <- vlite::resizeBed(
      bed,
      center = "region",
      upstream = ext,
      downstream = ext,
      genome= genome,
      ignore.strand = FALSE # symmetric extension -> no strand needed
    )
  }
  
  # Compute bins ----
  bins <- if(method=="nbins") {
    GenomicRanges::tile(GenomicRanges::GRanges(bed),
                        n = nbins)
  } else if(method=="slidingWindow") {
    GenomicRanges::slidingWindows(GenomicRanges::GRanges(bed),
                                  width= bins.width,
                                  step = steps.width)
  }
  bins <- as.data.table(bins)
  bins$group_name <- bins$width <- NULL
  
  # Expand list ----
  res <- data.table::copy(bed)
  res <- res[, setdiff(names(bed), names(bins)), with= F]
  res <- res[(bins$group)]
  bins$group <- NULL
  res <- cbind(bins, res)
  
  # Bin indices ----
  res[, bin.idx:= rowid(line.idx)]
  if(!ignore.strand)
    res[strand=="-", bin.idx:= rev(bin.idx), line.idx]
  
  # Order columns and return ----
  setcolorder(res, c("line.idx", "bin.idx"))
  return(res)
}
