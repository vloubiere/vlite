#' Motif enrichment analysis
#'
#' Compute motif enrichment between a set of regions and control regions
#'
#' @param counts data.table containing counts for the regions of interest.
#' @param control.counts data.table containing counts for control regions. Should have the same number of columns as 'counts' table.
#' @param names Convenient names to be used for aggregating and plotting. By default, returns motif_cluster matching motif_ID in vl_Dmel_motifs_DB_full.
#' @param plot Plot result?
#' @param padj.cutoff cutoff for plotting. Default= FALSE.
#' @param top.enrich Show only n top enriched motifs. Default= Inf (all).
#' @param min.counts The minimum number of counts required to call a hit. Default= 3L.
#' @param breaks Color breaks to be used. Defaults to range of filtered padj.
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj", "log2OR". Defaut= "padj".
#' @param add.motifs Should motif pwms be plotted?
#' @param cex.width Expansion factor for motif widths.
#' @param col Colors vector for bars.
#' @param xlab x label.
#' @param cex.height Expansion factor for motif heights.
#'
#' @examples
#' # Resize example peaks
#' SUHW <- resizeBed(vl_SUHW_top_peaks, genome = "dm3")
#' STARR <- resizeBed(vl_STARR_DSCP_top_peaks, genome = "dm3")
#'
#' # Generate same number of random regions
#' random <- controlRegionsBSgenome(bed= STARR, genome= "dm3")
#'
#' # Count JAPSPAR motifs (see below to use custom list of PWMs)
#' suhw <- vl_motifCounts(SUHW, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' starr <- vl_motifCounts(top_STARR, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' ctl <- vl_motifCounts(random, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#'
#' # Starting from sequence instead of bed file
#' seq <- getBSsequence(SUHW, genome= "dm3")
#' seq_suhw <- vl_motifCounts(seq, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' identical(suhw, seq_suhw)
#'
#' # Motifs can also be counted using a custom PWMatrixList, for example for promoter motifs:
#' prom_db <- readRDS("/groups/stark/almeida/data/motifs/CP_motifs/CP_motifs_PWM.rds")
#' prom <- prom_db$Pwms_log_odds
#' for(i in seq(prom))
#'   prom[[i]]@profileMatrix <- pwmPercToLog(prom_db$Pwms_perc[[i]]@profileMatrix)
#'
#' prom_motifs <- vl_motifCounts(STARR,
#'                                pwm_log_odds= prom,
#'                                genome= "dm3")
#'
#' # Compute enrichment at SUHW peaks, using random controls as background
#' pl <- vl_motifEnrich(suhw,
#'                       ctl,
#'                       plot= F)
#' pl[order(padj)][1:3] # Top motif is su(Hw)
#' # Plot
#' plot(pl,
#'      top.enrich= 3)
#'
#'
#' # Compute enrichment at STARR peaks and plot at the same time
#' enr <- vl_motifEnrich(starr,
#'                        ctl,
#'                        plot= T,
#'                        order= "log2OR",
#'                        padj.cutoff= 1e-5)
#' # Positive enrichments identify typical S2 enhancer motifs
#' plot(enr[log2OR>.5])
#'
#' @return DT of enrichment values which can be plot using ?plot.vl_GO_enr
#' @export
vl_motifEnrich <- function(counts,
                           control.counts,
                           names= vl_Dmel_motifs_DB_full[colnames(counts), motif_cluster, on= "motif_ID"],
                           plot= F,
                           padj.cutoff= 0.05,
                           top.enrich= Inf,
                           min.counts= 3L,
                           order= "padj",
                           breaks= NULL,
                           col= c("blue", "red"),
                           xlab = "Odd Ratio (log2)",
                           add.motifs= F,
                           cex.width= 1,
                           cex.height= 1)
{
  if(!is.data.table(counts))
    counts <- as.data.table(counts)
  if(!is.data.table(control.counts))
    counts <- as.data.table(control.counts)
  if(!all(sapply(counts, class)=="numeric"))
    stop("counts should only contain numeric values")
  if(!all(sapply(control.counts, class)=="numeric"))
    stop("control.counts should only contain numeric values")
  if(!is.null(names) && length(names)!=ncol(counts))
    stop("names should match ncol(counts)")

  # make obj ----
  obj <- rbindlist(list(set= counts,
                        control= control.counts),
                   idcol = T)

  # Melt ----
  obj <- melt(obj,
              id.vars = ".id",
              variable.name= "variable")

  # Test enrichment
  res <- obj[, {
    # Contingency table
    tab <- table(factor(.id=="set", levels = c(T, F)),
                 factor(value>0, levels = c(T, F)))
    # Check contingency table -> Fisher
    res <- fisher.test(tab)
    .(OR= res$estimate,
      pval= res$p.value,
      set_hit= sum(.id=="set" & value>0),
      set_total= sum(.id=="set"),
      ctl_hit= sum(.id=="control" & value>0),
      ctl_total= sum(.id=="control"))
  }, variable]

  # Add names
  res[, name:= names]
  res[is.na(name), name:= variable] # If names can't be found

  # padj...
  res[, padj:= p.adjust(pval, method = "fdr")]
  res[, log2OR:= log2(OR)]
  res$OR <- NULL

  # Order and save
  setcolorder(res, c("variable", "log2OR", "pval", "padj"))
  setattr(res, "class", c("vl_enr", "data.table", "data.frame"))
  if(plot)
  {
    DT <- plot.vl_enr(obj= res,
                      padj.cutoff= padj.cutoff,
                      top.enrich= top.enrich,
                      min.counts= min.counts,
                      order= order,
                      breaks= breaks,
                      xlab = xlab,
                      col = col)
    if(add.motifs)
      addMotifs(DT,
                cex.width= cex.width,
                cex.height= cex.height)
  }
  invisible(res)
}
