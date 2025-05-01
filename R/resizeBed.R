#' Resize Genomic Regions in a BED File
#'
#' @description
#' Adjusts the size of genomic coordinates.
#'
#' @param bed Input genomic ranges in a format compatible with ?importBed().
#' @param center Specifies the origin for resizing. Options are:
#' \itemize{
#'   \item `center`: Resize symmetrically around the midpoint of the region (default).
#'   \item `start`: Resize from the start of the region.
#'   \item `end`: Resize from the end of the region.
#'   \item `region`: Extend or shorten the entire region's boundaries.
#' }
#' @param upstream Integer. Number of bases to extend upstream from the specified origin. Default= 500.
#' @param downstream Integer. Number of bases to extend downstream from the specified origin. Default= 500L.
#' @param genome A BS genome name (e.g. "dm6", "mm10").
#' If specified, resized regions exceeding chromosome lengths will be clipped. Default= NULL.
#' @param ignore.strand Logical. If set to TRUE, upstream and downstream resizing is performed irrespective
#' of the strand. Default= FALSE.
#'
#' @details
#' Resized regions are returned as a gr data.table containing any additional columns from the input.
#'
#' @examples
#' # Import example BED file
#' bed <- importBed(c("chr2L:1000-2000:+", "chr2L:1000-2000:-"))
#'
#' # Resize using different centering options
#' resizeBed(bed, center = "start", upstream = 500, downstream = 1000)[]
#' resizeBed(bed, center = "center", upstream = 500, downstream = 1000)[]
#' resizeBed(bed, center = "end", upstream = 500, downstream = 1000)[]
#' resizeBed(bed, center = "region", upstream = 500, downstream = 1000)[]
#'
#' # Ignore strand during resizing
#' resizeBed(bed, center = "center", upstream = 500, downstream = 1000, ignore.strand = TRUE)[]
#'
#' @export
resizeBed <- function(bed,
                      center= "center",
                      upstream= 500L,
                      downstream= 500L,
                      genome= NULL,
                      ignore.strand= FALSE)
{
  # Import bed ----
  current <- importBed(bed)

  # Set strand ----
  if(ignore.strand | !"strand" %in% names(current))
    current[, strand:= "*"]

  # Get region start depending on strand ----
  if(center=="start"){
    current[, newStart:= ifelse(strand=="-", end, start)]
  } else if(center=="end") {
    current[, newStart:= ifelse(strand=="-", start, end)]
  } else if(center=="center") {
    current[, newStart:= floor(rowMeans(.SD)), .SDcols= c("start", "end")]
  } else if(center=="region") {
    current[, newStart:= start]
  }

  # Resize ----
  current[, c("newStart", "newEnd"):= {
    .(newStart-ifelse(strand=="-", downstream, upstream),
      newStart+ifelse(strand=="-", upstream, downstream))
  }]
  if(center=="region")
    current[, newEnd:= newEnd+(end-start)]

  # Modify bed and return ----
  res <- importBed(bed)
  res[, start:= current$newStart]
  res[, end:= current$newEnd]

  # Clip regions outside of chromosomes ----
  if(!is.null(genome)) {
    sizes <- getBSgenomeSize(genome = genome)
    res <- clipBed(res, sizes)
  }

  # Return ----
  return(res)
}
