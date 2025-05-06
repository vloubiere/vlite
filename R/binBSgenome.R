#' Bin BS genome
#'
#' Very fast using data.table
#'
#' @param genome BSgenome object to use. Example "dm3"
#' @param bins.width bins width default to 50bp
#' @param steps.width steps width separating each bin. default set to bins.width
#' @param restrict.seqnames If specified, bins are restricted to provided seqnames
#'
#' @examples
#' vl_binBSgenome(genome= "dm3",
#'                bins.width= 50)
#'
#' @return data.table containing bin coordinates
#' @export
binBSgenome <- function(genome,
                        bins.width= 50L,
                        steps.width= bins.width,
                        restrict.seqnames= NULL)
{
  # Retrive genome size
  gSize <- getBSgenomeSize(genome= genome)

  # Restrict to chromosomes
  if(!is.null(restrict.seqnames))
    gSize <- gSize[seqnames %in% restrict.seqnames]

  # Bin genome
  bins <- binBed(gSize,
                 bins.width = bins.width,
                 steps.width = steps.width,
                 ignore.strand = TRUE)

  # Return
  return(bins)
}
