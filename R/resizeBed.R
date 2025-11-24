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
  newStart <- if(center=="start"){
    current[, ifelse(strand=="-" & !ignore.strand, end, start)]
  } else if(center=="end") {
    current[, ifelse(strand=="-" & !ignore.strand, start, end)]
  } else if(center=="center") {
    current[, floor(rowMeans(.SD)), .SDcols= c("start", "end")]
  } else if(center=="region") {
    current[, start]
  } else
    stop("center should be one of 'center', 'start', 'end', 'region'.")

  # Compute width if relevant ----
  if(center=="region")
    width <- current[, end-start] # In this case, do not subtract 1!

  # Resize ----
  current$start <- newStart-ifelse(current$strand=="-" & !ignore.strand, downstream, upstream)
  current$end <- newStart+ifelse(current$strand=="-" & !ignore.strand, upstream, downstream)
  if(center=="region")
    current$end <- current$end+width

  # Clip regions outside of chromosomes ----
  if(!is.null(genome)) {
    sizes <- vlite::getBSgenomeSize(genome = genome)
    current <- vlite::clipBed(current, sizes)
  }

  # Return ----
  return(current)
}
