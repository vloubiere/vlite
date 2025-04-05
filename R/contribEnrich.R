#' Compute enrichment for each motif instance
#'
#' @param mot A list of motif created with !vl_motif_pos()
#' @param contrib.score A list of contribution scores, as in ?importContrib() output
#' @param seq.length The length of the tiles used in the model. Default= 1001L.
#' @examples
#' # Contrib oject ----
#' files <- list.files(paste0("db/", model, "/contributions/"),
#'                     ".h5$",
#'                     recursive = T,
#'                     full.names = T)
#' contrib <- importContrib(h5 = files)
#'
#' # Find motifs ----
#' mot <- vl_motif_pos(sequences = contrib,
#'                     pwm_log_odds = vl_motifs_DB_v2[collection=="flyfactorsurvey" & !is.na(Dmel), pwms_log_odds],
#'                     bg = "genome",
#'                     genome = "dm6",
#'                     p.cutoff = 1e-4,
#'                     collapse.overlapping = TRUE)
#'
#' # Motif enrichment ----
#' enr <- contribEnrich(mot = mot,
#'                      contrib.score = contrib$score,
#'                      seq.length = 1001L)
#'
#' @return A data.table with the first columns corresponding to sequence levels (factor) and each column to the matches for the corresponding motif.
#' @export
contribEnrich <- function(mot,
                          contrib.score,
                          seq.length= 1001L)
{
  # Remove unused levels if any ----
  mot[, seqlvls:= factor(seqlvls, unique(seqlvls))]
  if(length(levels(mot$seqlvls))!=length(contrib.score))
    stop("The number of levels in seqlevels in mot should matche the length(contrib.score)")

  # Scale score from 0 to 1 and compute means ----
  min.score <- min(unlist(contrib.score), na.rm = TRUE)
  range.score <- diff(range(unlist(contrib.score), na.rm = TRUE))
  scaled.score <- lapply(contrib.score, function(x) (x-min.score)/range.score)
  mean.scaled.score <- sapply(scaled.score, mean)

  # Retrieve regions for rdm controls ----
  regions <- mot[, .(seqnames= levels(seqlvls), start= 1, end= seq.length)]
  mot.widths <- mot[!is.na(mot.count), round(mean(rbindlist(ir)$width)), motif]
  mot[mot.widths, bins.width:= i.V1, on= "motif"]

  # Loop over motifs and compute enrichment ----
  res <- mot[!is.na(bins.width), {
    # Bin all sequences
    bins <- binBed(bed = regions,
                   bins.width = bins.width)
    # Select bins with correct width
    bins[, width:= end-start+1]
    bins <- bins[width==bins.width]
    bins[, seqlvls:= factor(seqnames, levels(seqlvls))]
    # Clean and add group
    bins$seqnames <- bins$binIDX <- NULL
    # For each motif
    .SD[, {
      # 2X random sampling
      sel <- sum(mot.count, na.rm= TRUE)*2
      set.seed(.GRP*sel)
      rdm <- bins[sample(.N, sel, replace = sel>.N)]
      # Motif instances
      mot <- .SD[!is.na(mot.count), ir[[1]], seqlvls]
      # Combine
      cmb <- list(motif= mot, rdm= rdm)
      cmb <- rbindlist(cmb, fill= TRUE, idcol = "group")
      # Extract motif scores
      cmb[, scores:= .(scaled.score[seqlvls])]
      cmb[, mot.scores:= .(.(scaled.score[[seqlvls]][start:end])), .(seqlvls, start, end)]
      # Compute enrichment
      cmb[, mean.score:= sapply(mot.scores, mean)]
      cmb[, log2OR:= log2(mean.score/mean.scaled.score[seqlvls])]
      # Compute OR and pval for each motif instance
      cmb[, pval:= {
        wilcox.test(mot.scores[[1]],
                    scores[[1]],
                    alternative= "greater")$p.value
      }, .(seqlvls, start, end)]
      cmb[, padj:= p.adjust(pval, method = "fdr")]
      # Fisher test
      .t <- table(motif= factor(cmb$group=="motif", c(TRUE, FALSE)),
                  signif= factor(cmb$padj<0.05, c(TRUE, FALSE)))
      OR.motif <- fisher.test(.t+1,
                              alternative = "greater")$estimate
      pval.motif <- fisher.test(.t,
                                alternative = "greater")$p.value
      # Select significant motifs coordinates
      sig.mot.coor <- cmb[group=="motif" & padj<0.05]
      sig.mot.coor <- sig.mot.coor[, .(seqlvls, start, end, width, mean.score, log2OR, padj)]
      # Return
      .(log2OR= log2(OR.motif),
        pval= pval.motif,
        sig.inst= .t[1,1],
        tot.inst= sum(.t[1,]),
        sig.rdm= .t[2,1],
        tot.rdm= sum(.t[2,]),
        sig.mot.coor= .(sig.mot.coor))
    }, motif]
  }, bins.width]
  # ADjusted p values
  res[, padj:= p.adjust(pval, "fdr")]

  # Clean ----
  clean <- res[, .(motif, log2OR, padj,
                   sig.inst, tot.inst, sig.rdm, tot.rdm,
                   sig.mot.coor)]

  # Add missing values ----
  all <- data.table(motif= unique(mot$motif))
  final <- merge(all, clean, by= "motif", all.x= TRUE, sort= FALSE)

  # Return ----
  return(final)
}
