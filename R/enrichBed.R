#' Assess overlap between peaks and regions
#'
#' @description
#' Test whether a set of regions significantly overlaps a set of peaks.
#'
#' @param regions Regions for which overlaps with peaks have to be assessed, in any format compatible with ?importBed.
#' @param peaks Peaks regions in any format compatible with ?importBed.
#' @param ctl.peaks Control peak regions with no overlaps with peaks, in any format compatible with ?importBed.
#' @param pseudocount The pseudocount to use to avoid 0s in fisher test estimates. Default= 0.5.
#' Note that the p-value will not be affected.
#' @param maxgap A single integer specifying the maximum gap allowed between 2 ranges for them to
#' be considered as overlapping. Default= -1.
#' @param minoverlap A single integer specifying the minimum overlap between 2 ranges for them to
#' be considered as overlapping. Default= 0.
#' @param ignore.strand If set to FALSE, only overlaps between regions that are on the same strand will be counted.
#' If set to TRUE (default), overlaps on both strands are counted.
#'
#' @return A list with elements:
#' \itemize{
#'   \item tab: contingency table.
#'   \item p.value: fisher test p.value.
#'   \item estimate: the fisher test odd ratio, based on the specified pseudocount.
#' }
#'
#' @examples
#' peaks <- data.table(seqnames= "chr2R",
#'                     start= seq(1000, 10000, 1000))
#' peaks[, end:= start+100]
#' ctl.peaks <- resizeBed(peaks, "end", -100, 200)
#' enrichBed(peaks = peaks[1:5],
#'           ctl.peaks = ctl.peaks,
#'           regions = peaks)
#'
#' @export
enrichBed <- function(regions,
                      peaks,
                      ctl.peaks,
                      pseudocount= 0.5,
                      maxgap= -1L,
                      minoverlap= 0L,
                      ignore.strand= TRUE)
{
  # Import ----
  peaks <- importBed(peaks)
  ctl.peaks <- importBed(ctl.peaks)
  regions <- importBed(regions)

  # Checks ----
  if(any(covBed(peaks, ctl.peaks)>0))
    stop("There should be no overlaps between peaks and ctl.peaks.")

  # Compute overlaps ----
  # peaks
  peaks.cov <- covBed(peaks, regions, ignore.strand = ignore.strand)>0
  peaks.ov <- data.table(peak= rep(T, nrow(peaks)), region= ifelse(peaks.cov, T, F))
  # ccontrol peaks
  ctl.peaks.cov <- covBed(ctl.peaks, regions, ignore.strand = ignore.strand)>0
  ctl.peaks.ov <- data.table(peak= rep(F, nrow(ctl.peaks)), region= ifelse(ctl.peaks.cov, T, F))

  # Compute enrichment ----
  # Contingency matrix
  enr <- rbind(peaks.ov, ctl.peaks.ov)
  mat <- table(peaks= factor(enr$peak, c(TRUE, FALSE)),
               regions= factor(enr$region, c(TRUE, FALSE)))
  # p.value
  .f <- fisher.test(mat, alternative = "greater")
  # Use pseudocount to avoid Inf
  if(any(mat==0)) {
    mat <- mat+log2OR.pseudocount
    .f$estimate <- (mat[1,1] * mat[2,2]) / (mat[2,1] * mat[1,2])
  }

  # Return ----
  return(list(tab= mat,
              p.value= .f$p.value,
              estimate= .f$estimate))
}
