#' Find Overlapping Regions Between Two Sets of Genomic Intervals
#'
#' @description
#' A wrapper around ?GenomicRanges::findOverlaps() that identifies, for each genomic range in a,
#' the overlapping regions in b.
#'
#' @param peaks Peaks regions in any format compatible with ?importBed.
#' @param ctl.peaks Control peak regions with no overlaps with peaks, in any format compatible with ?importBed.
#' @param regions Regions for which overlaps with peaks have to be assessed, in any format compatible with ?importBed.
#' @param pseudocount The pseudocount to use to estimate the enrichment using the ?fisher.test function. Default= 1.
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
#'   \item peaks.total: total number of peaks.
#'   \item peaks.overlaps: number of peaks overlapping regions.
#'   \item p.value: fisher test p.value.
#'   \item estimate: the fisher test odd ratio, based on the specified pseudocount.
#' }
#'
#' @examples
#'
#' @export
bedEnrich <- function(peaks,
                      ctl.peaks,
                      regions,
                      pseudocount= 1,
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
  enr <- rbind(peaks.ov, ctl.peaks.ov)
  tab <- table(factor(enr$peak, c(TRUE, FALSE)),
               factor(enr$region, c(TRUE, FALSE)))
  p.value <- fisher.test(tab, alternative = "greater")$p.value
  estimate <- fisher.test(tab+pseudocount, alternative = "greater")$estimate

  # Return ----
  return(list(peaks.total= nrow(peaks),
              peaks.overlaps= sum(peaks.cov),
              p.value= p.value,
              estimate= estimate))
}
