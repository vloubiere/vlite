#' Resize Genomic Regions in a BED File
#'
#' @description
#' Adjusts the size of genomic regions in a BED file by extending or contracting them
#' from a specified origin. Supports strand-aware resizing and flexible centering options.
#'
#' @param bed Input genomic ranges in a format compatible with `importBed()`. This can include:
#'   - A character vector of ranges (e.g., `"chr:start-end:strand"`)
#'   - A `GRanges` object
#'   - A `data.table` or `data.frame` with columns `seqnames`, `start`, `end`, and optionally `strand`
#'   - A file path to a `.bed` file
#' @param center Character. Specifies the origin for resizing. Options are:
#' \itemize{
#'   \item `"center"`: Resize symmetrically around the midpoint of the region (default).
#'   \item `"start"`: Resize from the start of the region.
#'   \item `"end"`: Resize from the end of the region.
#'   \item `"region"`: Extend or contract the entire region's boundaries.
#' }
#' @param upstream Integer. Number of bases to extend upstream from the specified origin. Default is `500L`.
#' @param downstream Integer. Number of bases to extend downstream from the specified origin. Default is `500L`.
#' @param ignore.strand Logical. If `TRUE`, the strand is ignored, and resizing always proceeds
#'   from the leftmost coordinate. Default is `FALSE`.
#'
#' @details
#' The function allows flexible resizing of genomic regions:
#' - **Centering options**:
#'   - `"center"`: Resizes symmetrically around the midpoint of the region.
#'   - `"start"`: Resizes relative to the start coordinate, respecting strand orientation.
#'   - `"end"`: Resizes relative to the end coordinate, respecting strand orientation.
#'   - `"region"`: Extends or contracts the entire region, maintaining its original boundaries.
#' - **Strand handling**:
#'   - When `ignore.strand = FALSE`, resizing respects the strand:
#'     * Positive strand: upstream is subtracted, downstream is added.
#'     * Negative strand: upstream is added, downstream is subtracted.
#'   - When `ignore.strand = TRUE`, resizing always proceeds from the leftmost coordinate.
#'
#' The function ensures that the resized regions are returned as a properly formatted `data.table`,
#' preserving any additional columns from the input.
#'
#' @return A `data.table` containing the resized genomic ranges with the following columns:
#' \itemize{
#'   \item `seqnames`: Chromosome or sequence name
#'   \item `start`: Resized start position
#'   \item `end`: Resized end position
#'   \item Additional columns from the input
#' }
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
    current[, newStart:= rowMeans(.SD), .SDcols= c("start", "end")]
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

  # Return ----
  return(res)
}
