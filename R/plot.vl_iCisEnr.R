#' Function method to plot motifEnrichRanks recovery curves
#'
#' @param enr An enrichment table as returned by ?vl_iCisEnrich.
#' @param NES.cutoff The NES cutoff to select enriched motifs. Default= 3.
#' @param simplify.clusters Should only the best motif per clsuter be returned? Default= FALSE.
#' @param main Title. Default= "Enriched motifs".
#' @param legend.title Default= "NES".
#' @param cluster.rows Default= TRUE.
#' @param cluster.cols Default= FALSE.
#' @param show.legend Default= TRUE.
#' @param legend.cex Default= 0.6.
#' @param show.row.clusters Default= FALSE.
#' @param show.numbers Default= TRUE/
#' @param numbers.cex Default= 0.4.
#' @param add.logos If a path to a .rds pwm db is provided, corresponding logos will be added
#' on the left of the heatmap. Default= "/zssd/scratch/vincent.loubiere/motifs_db/JASPAR_CORE/20260812_JASPAR2026_CORE_NON_REDUNDANT.rds".
#' @param logo.cex.width Width expansion factor for motif logos. Default= .7.
#' @param logo.cex.height Height expansion factor for motif logos. Default= .7.
#' @param ... Additional arguments will be passed to ?vl_heatmap.
#'
#' @returns A recovery plot
#' @export
#'
#' @examples
plot.vl_iCisEnr <- function(
    enr,
    NES.cutoff= 3,
    simplify.clusters= FALSE,
    main= "Enriched motifs",
    legend.title = "NES",
    cluster.rows = TRUE,
    cluster.cols = FALSE,
    show.legend= TRUE,
    legend.cex= 0.6,
    show.row.clusters= FALSE,
    show.numbers= TRUE,
    numbers.cex= 0.4,
    add.logos= "/zssd/scratch/vincent.loubiere/motifs_db/JASPAR_CORE/20260812_JASPAR2026_CORE_NON_REDUNDANT.rds",
    logo.cex.width= 0.7,
    logo.cex.height= 0.7,
    ...
) {
  # Filter enriched motifs
  sel <- as.data.table(enr[motif %in% motif[NES>=NES.cutoff]])
  # Simplify per cluster (max NES)
  if(simplify.clusters) {
    sel[is.na(cluster), cluster:= motif]
    sel <- sel[, .SD[motif==motif[which.max(NES)]], cluster]
  }
  sel[, motif:= factor(motif, unique(motif))]
  # Cast matrix
  NES <- dcast(sel, motif~group, value.var = "NES")
  adj.NES <- NES <- as.matrix(NES, 1)
  # If padjust is provided, adjust the NES for plotting (n.s. in white)
  if("padj" %in% names(sel)) {
    padj <- dcast(sel, motif~group, value.var = "padj")
    padj <- as.matrix(padj, 1)
    adj.NES[padj>0.05] <- 0
  }
  # Plot heatmap
  hm <- vl_heatmap(
    x = adj.NES,
    cluster.rows = if(isTRUE(cluster.rows)) NES else cluster.rows,
    cluster.cols = cluster.cols,
    show.row.clusters = show.row.clusters,
    show.numbers = if(isTRUE(show.numbers)) round(NES, 1) else FALSE,
    numbers.cex = numbers.cex,
    legend.title = legend.title,
    show.legend = show.legend,
    legend.cex = legend.cex,
    ...
  )
  # Add logos
  if(!isFALSE(add.logos)) {
    mot <- readRDS(add.logos)
    addSeqLogoEnrPlot(
      obj = hm$rows[, .(name, y= y.pos+.2)],
      ICM = mot$ICM[rownames(NES)],
      cex.width = logo.cex.width,
      cex.height = logo.cex.height
    )
  }
  
  # Return
  invisible(hm)
}