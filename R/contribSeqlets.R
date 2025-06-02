#' Return seqLets with high contribution scores
#'
#' @param contrib A contribution data.table imported with ?importContrib.
#' @param ext.size The size by which candidate peaks should be extended to retrieve
#' local background. Default= 50.
#' @param log2OR.cutoff Only the peaks with |log2OR| >= log2OR.cutoff are returned.
#' @param FDR.cutoff Only the peaks with FDR <= FDR.cutoff are returned.
#' @examples
#' dir <- "/groups/stark/vloubiere/projects/epiDeepCancer/db/model_PHD11//contributions/"
#' files <- list.files(dir, ".h5", recursive= T, full.names = T)
#' contrib <- vlite::importContrib(files)
#' peaks <- contribSeqlets(contrib)
#'
#' @return A gr data.table containing significant peaks, based on provided cutoffs.
#' @export
contribSeqlets <- function(contrib,
                           ext.size= 50,
                           log2OR.cutoff= 1,
                           FDR.cutoff= 0.05)
{
  # Checks ----
  if(nrow(contrib) != length(unique(contrib$name)))
    stop("contrib should contain unique 'name'(s).")

  # Melt contribution scores ----
  .m <- contrib[, .(score= score[[1]]), .(seqlvls= name)]
  .m[, start:= seq(.N), seqlvls]
  .m[, end:= start, seqlvls]

  # Make all scores positive (to compute OR) ----
  .m[, pos.score:= score-min(score, na.rm = T), seqlvls]

  # Identify candidate peaks using zscore ----
  .m[, zscore:= scale(score), seqlvls]
  .m[, peak:= !between(zscore, -1.5, 1.5)]
  .m[, peak_id:= .GRP, .(seqlvls, rleid(peak))]
  .m[, width:= start[.N]-start[1]+1, peak_id]
  # Cutoff on width
  .m[width<4, peak:= FALSE]

  # Split peak and background regions ----
  # Peaks
  peaks <- .m[(peak), .(
    start= min(start),
    end= max(end),
    mean.score= mean(score),
    pos.score= .(pos.score)
  ), .(seqlvls, peak_id)]
  # bg
  bg <- .m[(!peak), .(
    start= min(start),
    end= max(end),
    pos.score= .(pos.score)
  ), .(seqlvls, peak_id)]

  # Extend peaks to retrieve local background ----
  peaks[, ext.start:= start-ext.size]
  peaks[, ext.end:= end+ext.size]
  # Local bg
  peaks$ctl <- bg[peaks, .(.(unlist(pos.score))), .EACHI, on= c("seqlvls", "start<=ext.end", "end>=ext.start")]$V1
  if(any(lengths(peaks$ctl))==0)
    stop("Background regions could not be found for some peaks. Increase ext.size?")

  # Compute log2OR, pval, FDR ----
  peaks[, log2OR:= {
    mapply(function(x,y) {
      log2(mean(x)/mean(y))
    }, x= pos.score, y= ctl)
  }]
  peaks[, pval:= {
    mapply(function(x,y) {
      wilcox.test(x, y)$p.value
    }, x= pos.score, y= ctl)
  }]
  peaks[, FDR:= p.adjust(pval, method = "fdr")]

  # Filter and clean ----
  peaks <- peaks[FDR <= FDR.cutoff & abs(log2OR) >= log2OR.cutoff]
  res <- peaks[, .(seqlvls, start, end, mean.score, log2OR, FDR)]

  # Return ----
  return(res)
}
