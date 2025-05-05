#' Find Closest Genomic Features Between Two Sets of Regions
#'
#' @description
#' Identifies the closest genomic features in set B for each region in set A,
#' similar to `bedtools closest`. Supports strand-aware distance calculations,
#' multiple closest features, and distance filtering.
#'
#' @param a Query regions in any format compatible with ?importBed().
#' @param b Target regions in any format compatible with ?importBed().
#' @param n Number of closest features to report for each query region.
#' Default= 1L.
#' @param min.dist Minimum absolute distance allowed between features. Default= 0L,
#' meaning that both overlapping (dist = 0) and non-overlapping (dist > 0) features will be reported.
#' @param max.dist Maximum absolute distance allowed between features.
#'   Features farther than this are excluded. Default= Inf.
#' @param ignore.strand If set to FALSE, only closest features that are on the same strand will be considered.
#' If set to TRUE (default), closest features are reported regardless of their strand, and their relative
#' positions will be stored using the sign of the 'dist' column (see details).
#'
#' @details
#' **Distance Calculation**:
#' - For non-overlapping features:
#'   * Distance is the base pairs between feature boundaries
#'   * Negative distances indicate upstream features
#'   * Positive distances indicate downstream features
#' - For overlapping features:
#'   * Distance is 0
#'
#' **Strand Handling**:
#' - When `ignore.strand = FALSE`:
#'   * Only features on matching strands are considered
#'   * Distance sign follows strand orientation
#' - When `ignore.strand = TRUE`:
#'   * Strand is ignored for matching
#'   * Distance sign follows genomic coordinates
#'
#' **Distance Filtering**:
#' - `min.dist`: Excludes features closer than this distance
#' - `max.dist`: Excludes features farther than this distance
#' - Both filters use absolute distances (direction ignored)
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item idx.a: Integer index of region in set A
#'   \item idx.b: Integer index of closest region(s) in set B
#'   \item dist: Distance between matched regions:
#'     - Negative: B is upstream of A
#'     - Zero: Regions overlap
#'     - Positive: B is downstream of A
#' }
#'
#' @examples
#' # Create example regions
#' a <- importBed(c(
#'   "chr2L:10-20:+",    # Region A1
#'   "chr2L:10-20:-",    # Region A2
#'   "chr3L:10-20:+"     # Region A3
#' ))
#' b <- importBed(c(
#'   "chr2L:5-7:+",      # Region B1
#'   "chr2L:23-26:-",    # Region B2
#'   "chr2L:100-200:+",  # Region B3
#'   "chr2L:1000-2000:+" # Region B4
#' ))
#'
#' # Find single closest feature (strand-agnostic)
#' idx <- closestBed(a, b)
#' cbind(idx, a[idx$idx.a], b[idx$idx.b])
#'
#' # Find closest feature on same strand
#' idx <- closestBed(a, b, ignore.strand = FALSE)
#' cbind(idx, a[idx$idx.a], b[idx$idx.b])
#'
#' # Find three closest features
#' idx <- closestBed(a, b, n = 3)
#' cbind(idx, a[idx$idx.a], b[idx$idx.b])
#'
#' # Filter by distance
#' # Only features at least 10bp apart
#' idx <- closestBed(a, b, min.dist = 10)
#'
#' # Only features within 100bp
#' idx <- closestBed(a, b, max.dist = 100)
#'
#' # Features between 150-1000bp apart
#' idx <- closestBed(a, b, min.dist = 150, max.dist = 1000)
#'
#' @export
closestBed <- function(a,
                       b,
                       n= 1,
                       min.dist= 0L,
                       max.dist= Inf,
                       ignore.strand= TRUE)
{
  # Checks ----
  if(min.dist<0L | max.dist<0L)
    stop("min.dist and max.dist, should be positive integers")

  # Import bed files ----
  a <- importBed(a)
  b <- importBed(b)

  # Should strand be considered?
  .cols <- if(!ignore.strand && "strand" %in% names(a) & "strand" %in% names(b))
    c("seqnames", "strand") else
      "seqnames"

  # Compute closest ----
  idx <- b[a, {
    # Compute all distances
    dist <- fcase(x.start>i.end, as.integer(x.start-i.end),
                  x.end<i.start, as.integer(i.start-x.end),
                  is.na(x.start) | is.na(x.end), NA_integer_,
                  default= 0L)
    # Compute range depening on n
    range <- sort(unique(dist[between(dist, min.dist, max.dist)]))
    min <- range[1]
    max <- range[ifelse(length(range) > n, n, length(range))]
    # Select closest features
    sel <- between(dist, min, max)
    .(idx.a= .GRP,
      idx.b= .I[sel],
      dist= dist[sel])
  }, .EACHI, on= .cols]

  # Remove regions in a with no closest feature in b ----
  count <- uniqueN(idx[is.na(idx.b), idx.a])
  if(length(count))
    message(paste0(count, "/", nrow(a),
                   " region(s) in 'a' had no closest features in b with these parameters and were discarded"))
  idx <- na.omit(idx[, .(idx.a, idx.b, dist)])

  # Adjust distance of upstream features ----
  upstream <- a[idx$idx.a, start] > b[idx$idx.b, end]
  idx[upstream, dist:= -dist]
  if("strand" %in% names(a))
  {
    neg.strand <- a[idx$idx.a, strand]=="-"
    idx[neg.strand, dist:= -dist]
  }

  # Return ----
  return(idx)
}
