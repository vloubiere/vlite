#' Resize Genomic Regions
#'
#' @description
#' Resizes genomic regions by `upstream` and `downstream` nucleotides around a chosen origin. The final
#' width is `upstream + downstream + 1`, except when `center = "region"`, which extends the
#' original region boundaries. Coordinates (`start`/`end`) are 1-based inclusive as returned by `importBed()`.
#'
#' @param bed Input regions in any format compatible with ?importBed.
#' @param center Specifies the origin for resizing. Options are:
#' \itemize{
#'   \item `center`: Resize around the region midpoint. For even-length regions, chooses the left of the two
#'   central bases, whereas `GenomicRanges::resize(fix="center")` uses an interbase midpoint convention.
#'   \item `start`: Resize from the 5' end when ignore.strand= FALSE, otherwise genomic start.
#'   \item `end`: Resize from the 3' end when ignore.strand= FALSE, otherwise genomic end.
#'   \item `region`: Extend or contract the entire region's boundaries.
#' }
#' @param upstream Size of the extension upstream of the specified origin. Default= 500L.
#' @param downstream Size of the extension downstream of the specified origin. Default= 500L.
#' @param genome A BS genome name (e.g. "dm6", "mm10") used to clipped resized regions.
#' Default= NULL (no clipping).
#' @param ignore.strand If set to TRUE, resizing always proceeds from the leftmost (start) coordinates,
#' irrespective of the strand. If set to FALSE (default), upstream and downstream resizing respect
#' the strand (start= 5' and end= 3'); unstranded (`*`) features are treated as `+`.
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
  # Checks ----
  center <- match.arg(center, c("start", "end", "center", "region"), several.ok = FALSE)
  
  # Import bed for incapsulation ----
  current <- vlite::importBed(bed)
  
  # Extend whole region ----
  if(center=="region") {
    current[, start:= ifelse(ignore.strand | strand!="-", start-upstream, start-downstream)]
    current[, end:= ifelse(ignore.strand | strand!="-", end+downstream, end+upstream)]
  } else {
    # Get extended region start ----
    if(center=="start") {
      current[, start:= ifelse(ignore.strand | strand!="-", start-upstream, end-downstream)]
    } else if(center=="end") {
      current[, start:= ifelse(ignore.strand | strand!="-", end-upstream, start-downstream)]
    } else if(center=="center") {
      current[, start:= {
        as.integer((start + end) %/% 2  - ifelse(ignore.strand | strand!="-", upstream, downstream))
      }]
    }
    # Compute extended region end ----
    current[, end:= start+upstream+downstream]
  }
  
  # Clip regions outside of chromosomes ----
  if(!is.null(genome)) {
    sizes <- vlite::getBSgenomeSize(genome = genome)
    missing.chr <- unique(setdiff(current$seqnames, sizes$seqnames))
    if(length(missing.chr))
      warning(
        paste0(
          paste0(
            "The following chromosomes are missing from the genome file and will be removed (n=",
            sum(current$seqnames %in% missing.chr),
            "):\n"
          ),
          paste0(missing.chr, collapse = ", ")
        )
      )
    # if(any(current$seqnames %in% ))
    current <- suppressWarnings(vlite::clipBed(current, sizes)) # Warnings about start/end coordinates
  }
  
  # Sanity check ----
  if(any(current$start<1 | current[,start>end]))
    warning("Some regions with start<1 or start>end!")
  
  # Return ----
  return(current)
}
