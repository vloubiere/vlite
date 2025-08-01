#' Plot and Cluster Top Marker Genes as a Heatmap from Single-Cell Transcriptome Data
#'
#' Selects top marker genes per cluster based on multiple statistical cutoffs, then visualizes their expression as a heatmap.
#' Optionally clusters rows and/or orders columns, and allows highlighting of selected genes.
#'
#' @param dat A data.table containing marker gene statistics (see ?sc_computeMarkerGenes).
#' @param topN Number of top genes to select per cluster (default: Inf, i.e., all).
#' @param perc.cutoff Minimum percentage of expressing cells within the cluster (default: 10).
#' @param perc.diff.cutoff Minimum difference in percentage of expressing cells within vs. outside the cluster (default: 20).
#' @param order.var Variable to order genes by; one of "perc.diff", "auc", "logFC", or "log2OR" (default: "perc.diff").
#' @param padj.wilcox.cutoff Adjusted Wilcoxon p-value cutoff (default: 0.05).
#' @param logFC.cutoff Log fold change cutoff (default: 0).
#' @param padj.fisher.cutoff Adjusted Fisher p-value cutoff (default: 1).
#' @param log2OR.cutoff Log2 odds ratio cutoff (default: 1).
#' @param auc.cutoff AUC cutoff (default: 0.5).
#' @param alternative "greater" for upregulated markers, "smaller" for downregulated (default: "greater").
#' @param cluster.rows Logical; whether to cluster heatmap rows (default: TRUE).
#' @param cutree.rows Number of row clusters to cut the dendrogram into (default: 10).
#' @param order.cols Logical; whether to order columns by cluster and rank (default: FALSE).
#' @param show.col.clusters Should column cluster names be shown? Default= TRUE.
#' @param selection Character vector of gene symbols to highlight or force-include in the heatmap (default: character()).
#'
#' @return Invisibly returns a list containing the heatmap matrix and the table of top marker genes.
#' @export
sc_markerHeatmap <- function(
    dat,
    topN= 3,
    perc.cutoff= 10,
    perc.diff.cutoff= 20,
    order.var= "perc.diff",
    padj.wilcox.cutoff= 0.05,
    logFC.cutoff= 0.25,
    padj.fisher.cutoff= NULL,
    log2OR.cutoff= NULL,
    auc.cutoff= NULL,
    alternative= "greater",
    cluster.rows= TRUE,
    cutree.rows= NULL,
    order.cols= FALSE,
    show.col.clusters= TRUE,
    selection= character(),
    pdf.file= NULL
) {
  # Extract top markers ----
  top <- sc_topMarkers(
    dat= dat,
    perc.cutoff= perc.cutoff,
    perc.diff.cutoff= perc.diff.cutoff,
    padj.wilcox.cutoff= padj.wilcox.cutoff,
    logFC.cutoff= logFC.cutoff,
    padj.fisher.cutoff= padj.fisher.cutoff,
    log2OR.cutoff= log2OR.cutoff,
    auc.cutoff= auc.cutoff,
    order.var= order.var,
    topN= topN,
    selection = selection,
    alternative= alternative
  )

  # If marker genes were selected ----
  if(nrow(top)) {

    # Extract data ----
    mat <- dcast(dat[symbol %in% top$symbol],
                 cluster~symbol,
                 value.var = "perc",
                 drop= TRUE)
    mat <- as.matrix(mat, rownames= 1, drop= F)

    # Cluster rows ----
    if(isTRUE(cluster.rows)) {
      cl <- vl_heatmap(mat,
                       clustering.distance.rows = "spearman",
                       plot = F)
      mat <- mat[cl$rows[(order), name],]
    }

    # Order columns ----
    top[, cluster:= factor(cluster, rownames(mat))]
    if(order.cols) {
      setorderv(top, c("cluster", "rank"))
    } else if("selection" %in% names(top)) {
      setorderv(top, "selection")
      # Create a separate cluster at the end
      top[(selection), cluster:= "selection"]
    }
    mat <- mat[, unique(top$symbol), drop= F]

    # Columns annotations ----
    # annot <- dat[match(colnames(mat), symbol), Nat.cluster]
    # annot.col <- c("lightgrey", "chartreuse", "red", "gold", "royalblue1", "royalblue3", "blue")

    # Heatmap ----
    vl_heatmap(mat,
               cluster.rows = cluster.rows,
               clustering.distance.rows = "spearman",
               cutree.rows = cutree.rows,
               cluster.cols = top[, cluster[1], symbol, drop= F]$V1,
               # col.annotations = annot,
               # col.annotations.col = adjustcolor(annot.col, .5),
               # col.annotations.title = "Clusters Parreno",
               col= c("white", "red"),
               legend.title = "% exp. cells",
               show.col.clusters = show.col.clusters,
               col.gap.width = .5,
               row.gap.width = .5,
               legend.cex = .6,
               pdf.file = pdf.file)
  } else
    warning("No genes selected. Either relax cutoffs or provide a selection of genes.")
}
