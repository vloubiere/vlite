#' random region
#'
#' Sample random regions from BSgenome
#'
#' @param n Number of regions to sample
#' @param width Widths of regions to sample (length should either be 1 or equal n)
#' @param restrict.seqnames If specified, only the provided seqnames will be used
#'
#' @examples
#' randomRegionsBSgenome("dm3", 100, 1)
#'
#' @return data.table containing randomly sampled regions
#' @export
randomRegionsBSgenome <- function(genome,
                                  n,
                                  width= 1,
                                  restrict.seqnames= NULL)
{
  # Retrieve chromosome sizes
  gSize <- getBSgenomeSize(genome= genome)
  gSize[, width:= end-start+1]

  # Restrict seqnames
  if(!is.null(restrict.seqnames))
    gSize <- gSize[seqnames %in% restrict.seqnames]

  # Intitiate data.table
  rdm <- data.table(idx= seq(n),
                    width= width)

  # Random sampling
  rdm <- rdm[, {
    # Select chromosomes that are big enough
    .c <- gSize[.BY, on= "width>=width"]
    if(!nrow(.c))
      stop("Some widths are lager than all provided chromosomes!")
    # Adjust size to avoid going over chrom limits
    .c[, end:= end-width+1]
    # Randomly sample chromosomes based on their width
    idx <- sample(nrow(.c), size = .N, replace = TRUE, prob = .c$end)
    .c <- .c[(idx)]
    # Randomly sample start
    .c[, start:= sample(end, size = .N, replace = TRUE), .(seqnames, end)]
    # Compute end
    .c[, end:= start+width-1]
    # Return
    .c[, .(seqnames, start, end)]
  }, width]

  # Return result
  setcolorder(rdm, c("seqnames", "start", "end"))
  return(rdm)
}
