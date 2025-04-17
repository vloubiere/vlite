#' Collapse or Merge Overlapping Genomic Ranges
#'
#' @description
#' Merges genomic ranges that are within a specified distance of each other or overlap.
#' Provides options for strand-specific merging and index tracking of merged regions.
#'
#' @param bed Input genomic ranges in any format compatible with `importBed()`:
#' \itemize{
#'   \item Character vector of ranges ("chr:start-end[:strand]")
#'   \item GRanges object
#'   \item data.frame/data.table with required columns
#'   \item Path to a BED file
#' }
#' @param max.gap Integer. Maximum distance between features to be merged:
#' \itemize{
#'   \item Positive value: Merges features separated by ≤ max.gap bases
#'   \item Zero: Merges only touching or overlapping features
#'   \item Negative value: Requires overlap of ≥ |max.gap| bases to merge
#' }
#' Default is `0L` (merge touching or overlapping features).
#' @param return.idx.only Logical. If `TRUE`, returns indices of merged regions
#'   instead of collapsing them. Default is `FALSE`.
#' @param ignore.strand Logical. If `TRUE`, merges regions regardless of strand.
#'   If `FALSE`, only merges regions on the same strand. Default is `TRUE`.
#'
#' @details
#' The function performs iterative merging of genomic regions:
#'
#' **Merging Process**:
#' 1. Sorts regions by chromosome and coordinates
#' 2. Identifies adjacent or overlapping regions based on `max.gap`
#' 3. Iteratively merges regions until no more merging is possible
#'
#' **Strand Handling**:
#' - When `ignore.strand = TRUE`: Merges regions based only on genomic coordinates
#' - When `ignore.strand = FALSE`: Merges regions only if they are on the same strand
#'
#' **Output Options**:
#' - With `return.idx.only = FALSE`: Returns merged regions
#' - With `return.idx.only = TRUE`: Returns original regions with merge group indices
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
#'   "chr2R:1500-2500:-",  # Region 2 (overlaps Region 1)
#'   "chr2R:3000-4000:+"   # Region 3 (separate)
#' ))
#'
#' # Merge overlapping regions (strand-agnostic)
#' collapseBed(bed)
#'
#' # Merge regions considering strand
#' collapseBed(bed, ignore.strand = FALSE)
#'
#' # Get merge indices instead of merging
#' collapseBed(bed, return.idx.only = TRUE)
#'
#' # Merge regions within 500bp of each other
#' collapseBed(bed, max.gap = 500)
#'
#' # Merge only regions with at least 100bp overlap
#' collapseBed(bed, max.gap = -100)
#'
#' @export
collapseBed <- function(bed,
                        max.gap= 0L,
                        return.idx.only= F,
                        ignore.strand= TRUE)
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
