#' Randomly Sample Control Genomic Regions
#'
#' @description
#' Randomly sample genomic regions from a bed, matching specified widths and potentially avoiding overlaps
#' with other regions.
#'
#' @param bed Input genomic ranges, in any format compatible with ?importBed.
#' @param widths Integer vector specifying the width of the regions to sample. The length of the vector
#' determines the number of sampled regions.
#' @param no.overlaps Genomic ranges for which overlaps should be avoided dhring sampling, in any format compatible
#' with ?importBed. Default= NULL.
#' @param ignore.strand If genomic ranges are provided in no.overlaps, only overlaps on the same strand will be
#' avoided. If set to TRUE (default), overlaps on both strands are avoided.
#'
#' @return A data.table of sampled control regions.
#'
#' @examples
#' # Using widths
#' bed <- importBed(c("chr2L:1000-2000", "chr3R:1000-2000"))
#' randomRegionsBed(bed= bed,  widths= rep(10, 10))
#'
#' # Using a sample.bed file
#' no.overlaps <- importBed(c("chr2L:1800-2000", "chr3R:1800-2000"))
#' randomRegionsBed(bed= bed, widths= rep(10, 10), no.overlaps= no.overlaps)
#'
#' @export
randomRegionsBed <- function(bed,
                             widths= rep(100, 1000),
                             no.overlaps= NULL,
                             ignore.strand= TRUE)
{
  # Import bed ----
  bed <- importBed(bed)

  # Import sample.bed ----
  if(!is.null(no.overlaps)) {
    no.overlaps <- importBed(no.overlaps)
    bed <- subtractBed(a = bed,
                       b = no.overlaps,
                       ignore.strand = ignore.strand)
  }
  # Compute width after subtracting
  bed[, width:= end-start+1]

  # Random sampling ----
  widths <- data.table(width= widths)
  rdm <- widths[, {
    # Select sequences that are large enough
    sel <- bed[.BY, on= "width>=width", nomatch= NULL]
    if(nrow(sel)==0)
      stop("Some ranges are too big for the remaining regions.")
    # Resize based on sel width
    sel[, end:= end-width+1]
    # Randomly sample based on prob (see earlier)
    idx <- sample(x = seq(nrow(sel)),
                  size = .N,
                  replace = TRUE,
                  prob = sel[, end-start+1])
    .s <- sel[idx]
    # Compute new start and end
    .s[, new.start:= sample(x = start:end,
                            size = .N,
                            replace = TRUE), .(start, end)]
    .s[, new.end:= new.start+width-1]
    # Return
    .s[, .(seqnames, start= new.start, end= new.end, strand)]
  }, width]
  rdm$width <- NULL

  # Return result ----
  return(rdm)
}
