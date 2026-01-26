#' Find Closest Genomic Features Between Two Sets of Regions
#'
#' @description
#' A wrapper around `GenomicRanges::nearest()` and `GenomicRanges::nearestKNeighbors()` that computes,
#' for each genomic range in a, the nearest genomic ranges in b and the stranded distance between them.
#'
#' @param a Query regions in any format compatible with ?importBed.
#' @param b Target regions in any format compatible with ?importBed.
#' @param k If specified, the k nearest neighbors are computed. If set to NA (default), only the closest
#' regions are returned, including ties.
#' @param ignore.strand If set to FALSE, only closest features that are on the same strand are reported.
#' If set to TRUE (default), closest features are reported on both strands. See details.
#'
#' @details
#' **Distance Calculation**:
#' - Features that have no nearest neighbor (e.g., chromosome missing in b) are not returned
#' - For overlapping features, a distance of 0 is returned.
#' - For non-overlapping features, the genomic distance between their closest boundaries are returned:
#'   * Negative distances indicate upstream features.
#'   * Positive distances indicate downstream features.
#'   * Unstranded (*) features are treated as `+`.
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item idx.a: line index of the regions in a.
#'   \item idx.b: line index of the closest region(s) in b.
#'   \item dist: distance between closest regions (see details).
#' }
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr2R:100-200:+", "chr2R:100-200:-"))
#' b <- importBed(c("chr2R:300-400:+", "chr2R:300-400:-", "chr2R:3000-3100:+"))
#'
#' # Find single closest features
#' closestBed(a, b)[]
#'
#' Only consider features that are on the same strand
#' closestBed(a, b, ignore.strand = F)[]
#'
#' # Return all features at the second closest distance (including ties):
#' closestBed(a, b, k= 3)[]
#' closestBed(a, b, k= 3, ignore.strand= FALSE)[]
#'
#' @export
closestBed <- function(a,
                       b,
                       k= NA,
                       ignore.strand= TRUE)
{
  # Checks ----
  if(!is.na(k) && (k<1L | k %% 1 != 0))
    stop("If specified, k should be a positive integer")

  # Import a and b for encapsulation ----
  a <- vlite::importBed(a)
  b <- vlite::importBed(b)

  # Find nearest regions (n=1) ----
  if(is.na(k)) {
    closest <- suppressWarnings(
      GenomicRanges::nearest(
        GenomicRanges::GRanges(a),
        GenomicRanges::GRanges(b),
        select= "all",
        ignore.strand= ignore.strand
      )
    )
    closest <- as.data.table(closest)
    setnames(closest, c("idx.a", "idx.b"))
  } else {

    # Find k nearest neighbors ----
    a[, idx.a:= .I]
    b[, idx.b:= .I]
    # For each chromosome shared between a and b
    common.chr <- intersect(unique(a$seqnames), unique(b$seqnames))
    closest <- list()
    for(chr in common.chr) {
      .a <- a[seqnames==chr]
      .b <- b[seqnames==chr]
      # Find k nearest neighbors
      .cl <- GenomicRanges::nearestKNeighbors(
        GRanges(.a),
        GRanges(.b),
        k= k,
        ignore.strand= ignore.strand
      )
      .cl <- data.table(
        idx.a= rep(seq(.cl), lengths(.cl)),
        idx.b= unlist(.cl)
      )
      .cl[, idx.a:= .a$idx.a[idx.a]]
      .cl[, idx.b:= .b$idx.b[idx.b]]
      closest[[chr]] <- .cl
    }
    # Bind all chromosomes
    closest <- rbindlist(closest)
  }

  # Compute distances ----
  closest[, dist:= fcase(
    b$end[idx.b]<a$start[idx.a], b$end[idx.b]-a$start[idx.a],
    a$end[idx.a]<b$start[idx.b], b$start[idx.b]-a$end[idx.a],
    default = 0L
  )]

  # Upstream feature -> neg distance ----
  closest[, dist:= ifelse(a$strand[idx.a]=="-", -dist, dist)]

  # Return ----
  return(closest)
}
