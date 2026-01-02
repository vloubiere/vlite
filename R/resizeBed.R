#' Resize Genomic Regions in a BED File
#'
#' @description
#' Resize genomic ranges.
#'
#' @param bed Input regions in any format compatible with ?importBed.
#' @param center Specifies the origin for resizing. Options are:
#' \itemize{
#'   \item `center`: Resize symmetrically around the midpoint of the region (default).
#'   \item `start`: Resize from the start of the region.
#'   \item `end`: Resize from the end of the region.
#'   \item `region`: Extend or contract the entire region's boundaries.
#' }
#' @param upstream Size of the extension upstream of the specified origin. Default= 500L.
#' @param downstream Size of the extension downstream of the specified origin. Default= 500L.
#' @param genome A BS genome name (e.g. "dm6", "mm10") used to clipped resized regions.
#' Default= NULL (no clipping).
#' @param ignore.strand If set to TRUE, resizing always proceeds from the leftmost coordinates,
#' irrespective of the strand. If set to FALSE (default), upstream and downstream resizing respect
#' the strand.
#'
#' @return
#' A data.table containing resized regions.
#'
#' @examples
#' # Import example BED file
#' bed <- importBed(c("chr2L:1000-2000:+", "chr2L:1000-2000:-"))
#'
#' # Resize using different centering options
#' resizeBed(bed, center = "start", upstream = 100, downstream = 200)[]
#' resizeBed(bed, center = "center", upstream = 100, downstream = 200)[]
#' resizeBed(bed, center = "end", upstream = 100, downstream = 200)[]
#' resizeBed(bed, center = "region", upstream = 100, downstream = 200)[]
#'
#' # Ignore strand during resizing
#' resizeBed(bed, center = "center", upstream = 100, downstream = 200, ignore.strand = TRUE)[]
#'
#' @export
resizeBed <- function(bed,
                      center= "center",
                      upstream= 500L,
                      downstream= 500L,
                      genome= NULL,
                      ignore.strand= FALSE)
{
  # Import bed for incapsulation ----
  current <- vlite::importBed(bed)

  # Get region anchor depending on strand ----
  if(center=="start") {
    current[, start:= ifelse(!ignore.strand && strand=="-", end-(downstream+1), start-upstream)]
  } else if(center=="end") {
    current[, start:= ifelse(!ignore.strand && strand=="-", start-downstream, end-(upstream+1))]
  } else if(center=="center") {
    current[, start:= {
      ifelse(!ignore.strand && strand=="-", floor(rowMeans(.SD))-(downstream+1), floor(rowMeans(.SD))-upstream)
    }, .SDcols= c("start", "end")]
  } else if(center=="region") {
    current[, start:= ifelse(!ignore.strand && strand=="-", end+upstream, start-downstream)]
  } else
    stop("center should be one of 'center', 'start', 'end', 'region'.")
  current[, end:= start+upstream+downstream]

  # Clip regions outside of chromosomes ----
  if(!is.null(genome)) {
    sizes <- vlite::getBSgenomeSize(genome = genome)
    current <- vlite::clipBed(current, sizes)
  }

  # Sanity check ----
  if(any(current$start<1 | current[,start>end]))
    warning("Some regions with start<1 or start>end!")
  
  # Return ----
  return(current)
}
