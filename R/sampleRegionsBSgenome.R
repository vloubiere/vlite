#' Generate control regions
#'
#' Generate regions with similar dispersion, seqnames and widths distribution than the provided bed file
#'
#' @param bed Bed file used to produce similar control
#' @param genome BSgneome object to use. ex: "dm3", "dm6"
#' @param no.overlap If set to TRUE, avoids overlap between control sequences and the original bed file. Default= FALSE.
#'
#' @examples
#' controlRegionsBSgenome(vl_SUHW_top_peaks, "dm3")
#'
#' @return data.table containing control regions
#' @export
controlRegionsBSgenome <- function(bed,
                                   genome,
                                   no.overlap= FALSE)
{
  # Check bed input
  bed <- importBed(bed)
  bed[, width:= end-start+1]

  # Retrieve chromosome sizes
  gSize <- getBSgenomeSize(genome= genome)
  gSize[, width:= end-start+1]

  # Check that all bed chromosomes are present
  if(any(!bed$seqnames %in% gSize$seqnames))
    stop("Some of the seqnames within bed regions could not be found in the BSgenome.")

  # Subtract overlapping regions
  if(no.overlap)
    gSize <- subtractBed(gSize, bed)

  # Random sampling
  bed[, start:= {
    # Select contigs from the same chromosome
    .c <- gSize[.BY, on= c("seqnames", "width>=width")]
    if(!nrow(.c))
      stop("Some regions are larger than non-overlapping contigs!")
    # Adjust to avoid going over contigs limits
    .c[, end:= end-width+1]
    # Compute new width
    .c[, width:= end-start+1]
    # Randomly sample contigs based on their width
    idx <- sample(nrow(.c), .N, replace = TRUE, prob = .c$width)
    .c <- .c[idx]
    # Randomly sample start
    .c[, sample(start:end, .N, replace = TRUE), .(start, end)]$V1
  }, .(seqnames, width)]
  bed[, end:= start+width-1]

  # Return
  sel <- intersect(c("seqnames", "start", "end", "strand"),
                   names(bed))
  return(bed[, ..sel])
}
