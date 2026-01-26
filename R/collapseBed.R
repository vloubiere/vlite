#' Collapse or Merge Overlapping Genomic Ranges
#'
#' @description
#' A wrapper around `GenomicRanges::reduce()` that merges overlapping genomic ranges.
#'
#' @param bed Input regions in any format compatible with ?importBed.
#' @param min.gapwidth Minimum distance between features to be merged. Default= 1L.
#' \itemize{
#'   \item Positive value: Merges features separated by ≤ min.gapwidth bases.
#'   \item Zero: touching yet non-overlapping features will not be merged.
#' }
#' @param return.idx.only If set to TRUE, returns the indices of merged regions
#'   instead of collapsing them. Default= FALSE.
#' @param ignore.strand If set to FALSE, only overlapping regions that are on the same strand are merged.
#' If set to TRUE (default), regions are merged irrespective of their strand, which will be set to `*`.
#'
#' @return
#' If return.idx.only = FALSE: a gr data.table with columns:
#' \itemize{
#'   \item seqnames: chromosome or sequence name.
#'   \item start: start position of merged region.
#'   \item end: end position of merged region.
#'   \item strand: strand (*= unstranded).
#' }
#' If return.idx.only = TRUE: a run-length type id indicating overlapping regions belong
#' that can be merged.
#'
#' @examples
#' # Create example regions
#' bed <- importBed(c("chr2R:1000-2000:+", "chr2R:2001-2500:+", "chr2R:2300-3000:-", "chr2R:3500-4000:+"))
#'
#' # Merge overlapping or touching regions
#' collapseBed(bed)
#'
#' # Only merge overlapping regions
#' collapseBed(bed,  min.gapwidth= 0)
#'
#' # Only merge regions that are on the same strand
#' collapseBed(bed, ignore.strand = FALSE)
#'
#' # Get merge indices instead of merging
#' collapseBed(bed, ignore.strand = FALSE, return.idx.only = TRUE)
#'
#' # Merge regions within 500bp of each other
#' collapseBed(bed, min.gapwidth = 500)
#' collapseBed(bed, min.gapwidth = 500, ignore.strand = FALSE)
#' collapseBed(bed, min.gapwidth = 500, ignore.strand = FALSE, return.idx.only = TRUE)
#'
#' @export
collapseBed <- function(bed,
                        min.gapwidth= 1L,
                        return.idx.only= FALSE,
                        ignore.strand= TRUE)
{
  # Import for incapsulation ----
  bed <- vlite::importBed(bed)
  
  # Reduce ----
  gr <- GenomicRanges::GRanges(bed)
  coll <- GenomicRanges::reduce(
    gr,
    ignore.strand= ignore.strand,
    min.gapwidth = min.gapwidth
  )
  
  if(return.idx.only)
  {
    # Return run-length id ----
    idx <- GenomicRanges::findOverlaps(
      gr,
      coll,
      ignore.strand= ignore.strand
    )
    idx <- subjectHits(idx)
    return(idx)
  } else {
    
    # Order and return collapsed ranges ----
    res <- vlite::importBed(coll)
    setorderv(res,
              c("seqnames", "start", "end"))
    return(res)
  }
}
