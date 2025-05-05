#' Collapse or Merge Overlapping Genomic Ranges
#'
#' @description
#' Merges overlapping genomic ranges.
#'
#' @param bed Input genomic ranges in any format compatible with ?importBed().
#' @param max.gap Integer. Maximum distance between features to be merged:
#' \itemize{
#'   \item Positive value: Merges features separated by ≤ max.gap bases
#'   \item Zero: Merges only touching or overlapping features
#'   \item Negative value: Requires overlap of ≥ |max.gap| bases to merge
#' }
#' Default= 0L (merge touching or overlapping features).
#' @param return.idx.only If set to TRUE, returns the indices of merged regions
#'   instead of collapsing them. Default= FALSE.
#' @param ignore.strand If set to FALSE, only overlapping regions that are on the same strand
#' are merged. If set to TRUE (default), regions are merged regardless of their strand.
#'
#' @return Depends on `return.idx.only`:
#' \itemize{
#'   \item If `FALSE`: A data.table of merged regions with columns:
#'     - seqnames: Chromosome name
#'     - start: Start position of merged region
#'     - end: End position of merged region
#'     - strand: Strand (if present in input and `ignore.strand = FALSE`)
#'   \item If `TRUE`: A vector of indices indicating which regions were merged together
#' }
#' 
#' @examples
#' # Create example regions
#' bed <- importBed(c(
#'   "chr2R:1000-2000:+",  # Region 1
#'   "chr2R:1500-2500:+",  # Region 2 (overlaps Region 1)
#'   "chr2R:1500-2500:-",  # Region 3 (overlaps Region 1, but different strand)
#'   "chr2R:3000-4000:+"   # Region 4 (500 bp away from region 2)
#' ))
#'
#' # Merge overlapping regions with the same strand (3 merge regions)
#' collapseBed(bed)
#' 
#' # Get merge indices instead of merging ()
#' collapseBed(bed, return.idx.only = TRUE)
#' 
#' # Ignore strand (2 merged regions)
#' collapseBed(bed, ignore.strand = TRUE)
#'
#' # Merge regions within 450/500bp of each other
#' collapseBed(bed, max.gap = 450) # 3 merged regions
#' collapseBed(bed, max.gap = 500) # 2 merged regions
#'
#' # Merge only regions with at least 600/100bp overlap
#' collapseBed(bed, max.gap = -600) # 4 merged regions
#' collapseBed(bed, max.gap = -100) # 3 merged regions
#'
#' @export
collapseBed <- function(bed,
                        max.gap= 0L,
                        return.idx.only= FALSE,
                        ignore.strand= FALSE)
{
  # Import Bed ----
  bed <- importBed(bed)

  # Define groups that can be collapsed ----
  group_cols <- if(!ignore.strand && "strand" %in% names(bed))
    c("seqnames", "strand") else
      "seqnames"

  # Create copy ----
  ordBed <- data.table::copy(bed)
  ordBed[, bed.idx:= .I]
  ordBed[, idx:= .I]

  # Collapse ----
  setorderv(ordBed,
            c("seqnames", "start", "end"))
  before <- nrow(ordBed)
  after <- 0L
  while(after<before)
  {
    before <- nrow(ordBed)
    ordBed[, idx:= cumsum(!(start-max.gap-1L <= c(start[1], end[-(.N)]))), group_cols]
    # Make index >0 and avoid duplicated values on + and - strands
    ordBed[, idx:= .GRP, c(group_cols, "idx")]
    ordBed <- ordBed[, .(
      start= min(start),
      end= max(end),
      bed.idx= .(unique(unlist(bed.idx))) # original bed line indexes
    ), c(group_cols, "idx")]
    after <- nrow(ordBed)
  }

  # Move strand column to the end ----
  if("strand" %in% group_cols)
    setcolorder(ordBed, "strand", after= "end")

  if(return.idx.only)
  {
    # Return index ----
    idx <- rep(ordBed$idx, lengths(ordBed$bed.idx))[unlist(ordBed$bed.idx)]
    return(idx)
  }else
  {
    # Or collapsed table ----
    return(ordBed[, !c("idx", "bed.idx")])
  }
}
