#' Clip Genomic Regions to Defined Limits
#'
#' @description
#' Clips genomic regions in a to the boundaries defined by regions in b.
#' This function is used for restricting regions in a to valid genomic
#' boundaries, such as chromosome limits or other predefined regions.
#'
#' @param a Query regions in any format compatible with ?importBed().
#' @param b Target regions in any format compatible with ?importBed().
#' @param min.width Integer specifying the minimum width required for clipped regions.
#' Regions smaller than this are discarded. Default= 1L.
#' @param ignore.strand If set to FALSE, only clips at boundaries that are on the same strand.
#' If set to TRUE (default), clips boundaries regardless of their strands.
#'
#' @return A gr data.table containing the remaining portions of a, after clipping them using the boundaries
#' defined in b.
#'
#' @examples
#' # Create example regions
#' a <- importBed(c("chr1:100-500:+"))
#' b <- importBed(c("chr1:200-300:+", "chr1:400-450:-", "chr1:425-475:-"))
#'
#' # Basic example
#' clipBed(a, b)
#'
#' # Require minimum width of 100bp
#' clipBed(a, b, min.width= 100)
#' 
#' # Strand-specific
#' clipBed(a, b, ignore.strand = FALSE)
#'
#' @export
clipBed <- function(a,
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
  clip <- overlapBed(a,
                     coll,
                     ignore.strand = ignore.strand,
                     all.a = FALSE)

  # Overlaps coordinates correspond to clipped coor ----
  setnames(clip,
           c("overlap.start", "overlap.end"),
           c("start", "end"))

  # Select clipped regions larger than min.width ----
  clip <- clip[end-start+1 >= min.width]

  # cbind `a` ----
  clip <- cbind(a[clip$idx.a, .(seqnames)],
                clip[, .(start, end)],
                a[clip$idx.a, !c("seqnames", "start", "end")])

  # Return ----
  return(clip)
}
