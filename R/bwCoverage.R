#' bw coverage
#'
#' Compute mean signal from a bigwig file at a set of coordinates. Regions with no overlap will return NA.
#'
#' @param bed BED file input compatible with '?importBed', representing the genomic ranges for which mean bw signal will be returned.
#' @param bw Path to a target bw file.
#'
#' @examples
#' # Example track
#' bw <- "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/ATAC_PH18_merge.bw"
#'
#' Compute coverage
#' bed <- importBed(c("chr2R:10000000-10010000", "chr3R:10000000-10005000", "nonExistingChr:1000-2000"))
#' bwCoverage(bed= bed, bw= bw)
#'
#' @export
bwCoverage <- function(bed,
                       bw)
{
  # Checks ----
  if(length(bw) != 1)
    stop("length(bw) != 1")
  if(!file.exists(bw))
    stop("bw file does not exist!")

  # Hard copy for incapsulation, select cols and compute width ----
  regions <- importBed(bed)
  regions <- regions[, .(seqnames, start, end, width= end-start+1)]
  # Transform to Granges
  gr <- GenomicRanges::GRanges(na.omit(regions))

  # Import bw ----
  sel <- rtracklayer::BigWigSelection(gr, "score")
  var <- rtracklayer::import.bw(bw, selection= sel)
  var <- data.table::as.data.table(var)

  # Overlap ----
  res <- overlapBed(a= regions,
                    b= var,
                    ignore.strand = TRUE, # .bw files are not stranded
                    all.a = TRUE)
  # Retrieve regions widths (file a)
  res[, width:= regions[idx.a, width]]
  # Retrieve bins scores (file b)
  res[, score:= var[idx.b, score]]

  # Compute mean signal. *
  meanSig <- res[, sum(score*overlap.width)/width, keyby= .(idx.a, width)]$V1
  # * .bw files do not store NAs -> only regions with no overlaps will return NA

  # Return
  return(meanSig)
}
