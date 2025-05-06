#' Subtract Genomic Regions from Reference Intervals
#'
#' @description
#' Resize regions in a by subtracting regions in b in order to remove overlaps between
#' the two. If regions in a fully overlap region(s) in b, they will therefore be removed.
#'
#' @param a Query regions in any format compatible with ?importBed().
#' @param b Target regions in any format compatible with ?importBed().
#' @param min.width Integer specifying the minimum width required for resized regions.
#' Regions smaller than this are discarded. Default= 1L.
#' @param ignore.strand If set to FALSE, only subtracts overlapping features that are on the same strand.
#' If set to TRUE (default), subtracts overlapping feature(s) regardless of their strand(s).
#'
#' @return A gr data.table containing the remaining portions of a, after subtracting the regions
#' defined in b.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-500:-"))
#' b <- importBed(c("chr1:200-300:+", "chr1:400-450:-", "chr1:425-475:-"))
#'
#' # Basic example
#' subtractBed(a, b)
#' 
#' # Minimum width of 50bp
#' subtractBed(a, b, min.width = 50L)
#' 
#' # Only subtract regions with similar strand
#' subtractBed(a, b, ignore.strand= FALSE)
#'
#' @export
subtractBed <- function(a,
                        b,
                        min.width= 1L,
                        ignore.strand= TRUE)
{
  # Hard copy for incapsulation ----
  a <- importBed(a)

  # Collapse features to subtract ----
  coll <- collapseBed(b,
                      ignore.strand = ignore.strand)

  # Compute overlaps ----
  ov <- overlapBed(a,
                   coll,
                   ignore.strand = ignore.strand,
                   all.a = TRUE)
  # Subtract ----
  sub <- ov[, {
    .(seqnames= a$seqnames[idx.a],
      start= na.omit(c(a$start[idx.a], overlap.end+1)),
      end= na.omit(c(overlap.start-1, a$end[idx.a])))
  }, idx.a]

  # Select subtracted regions larger than min.width ----
  sub <- sub[end-start+1 >= min.width]

  # cbind `a` and remove index ----
  sub <- cbind(sub[, !"idx.a"],
               a[sub$idx.a, !c("seqnames", "start", "end")])

  # Return ----
  return(sub)
}
