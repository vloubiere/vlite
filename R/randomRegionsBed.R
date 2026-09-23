#' Randomly Sample Control Genomic Regions
#'
#' @description
#' Randomly sample genomic regions from a bed with replacement, matching specified widths and potentially
#' avoiding overlaps with another set or regions.
#'
#' @param bed Input genomic ranges (in any format compatible with ?importBed) from which random regions
#' are drawn.
#' @param widths Integer vector specifying the width of the regions to sample. The length of the vector
#' determines the number of sampled regions.
#' @param no.overlaps Genomic ranges for which overlaps should be avoided during sampling, in any format compatible
#' with ?importBed. Default= NULL.
#' @param ignore.strand If genomic ranges are provided in no.overlaps, only overlaps on the same strand will be
#' avoided. If set to TRUE (default), overlaps on both strands are avoided.
#'
#' @return A data.table of sampled control regions.
#'
#' @examples
#' # Draw regions from dm6 canonical chromosomes
#' dm6.sizes <- getBSgenomeSize(genome= "dm6")
#' dm6.sizes <- dm6.sizes[seqnames %in% c("chr2L", "chr2R")]
#' randomRegionsBed(bed= dm6.sizes,  widths= rep(10, 10))
#'
#' # Avoid overlaps with specific region
#' randomRegionsBed(bed= dm6.sizes, widths= rep(10, 10), no.overlaps= dm6.sizes[2])
#'
#' @export
randomRegionsBed <- function(bed,
                             widths= rep(100, 1000),
                             replace= T,
                             no.overlaps= NULL,
                             ignore.strand= TRUE)
{
  # Import bed ----
  bed <- importBed(bed)
  
  # Subtract regions specified in no.overlaps ----
  if(!is.null(no.overlaps)) {
    no.overlaps <- importBed(no.overlaps)
    bed <- subtractBed(
      a = bed,
      b = no.overlaps,
      ignore.strand = ignore.strand
    )
  }
  
  # Compute prob. based on width (after subtraction) ----
  bed[, prob:= end-start+1]
  # Sanity check
  if(any(widths>max(bed$prob)))
    stop("Some widths are wider than any remaining region.")
  if("width" %in% names(bed))
    bed$width <- NULL
  
  # Sample input regions that are wide enough for sampling ----
  rdm <- data.table(width= widths)
  rdm <- rdm[, {
    # Select based on size
    sel <- bed$prob>=width
    # Sample
    idx <- sample(
      x = which(sel),
      prob = bed$prob[sel],
      size= .N,
      replace = replace
    )
    bed[idx]
  }, width]
  
  # Randomly sample starting coordinates ----
  # Compute number of possible starts (increments)
  rdm[, npos:= prob-width+1]
  # Sample
  rdm[, add:= sample.int(npos, size= .N, replace= T)-1, npos]
  rdm[, new.start:= start+add]
  rdm <- rdm[, .(seqnames, start= new.start, end= new.start+width-1, strand, width)]
  
  # Return result ----
  return(rdm)
}
