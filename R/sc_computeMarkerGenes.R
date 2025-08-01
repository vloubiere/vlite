#' Compute Marker Genes for Selected Clusters
#'
#' Identifies marker genes for specified clusters in a Seurat object using Wilcoxon rank-sum test (via presto) and, optionally, Fisher's exact test.
#'
#' @param dat A Seurat object (must contain an SCT assay).
#' @param select.clusters Character vector of cluster names to analyze. If missing, all clusters are used.
#' @param add.fisher Logical; if TRUE, adds Fisher's exact test statistics for each gene/cluster.
#'
#' @return A data.table with marker gene statistics for each cluster, including average expression, log fold change, adjusted p-values, AUC, and
#' optionally Fisher's test results.
#'
#' @details
#' - Uses the SCT assay's data slot for expression values.
#' - Only clusters in `select.clusters` are analyzed if provided.
#' - Fisher's test is performed on active/inactive cell counts if `add.fisher = TRUE`.
#'
#' @import presto data.table
#' @export
sc_computeMarkerGenes <- function(
    dat,
    select.clusters= NULL,
    add.fisher= F
    )
{
  # Checks
  if(class(dat)!="Seurat")
    stop("dat should be a seurat object.")
  # Import data ----
  cl <- Idents(dat)
  SCT <- dat@assays$SCT$data

  # Subset selected clusters ----
  if(!is.null(select.clusters)) {
    select.clusters <- unique(select.clusters)
    if(!all(select.clusters %in% levels(cl)))
      stop("select.clusters missing from levels(Idents(dat))")
    # Simplify clusters
    cl[!cl %in% select.clusters] <- NA
    cl <- factor(cl, make.unique(c(select.clusters, "non.selected")))
    cl[is.na(cl)] <- rev(levels(cl))[1]
  }

  # Compute avgExpr, wilcox and auc ----
  stats <- presto::wilcoxauc(SCT, y= cl)
  stats <- as.data.table(stats)
  stats$statistic <- stats$pval <- NULL
  setnames(stats,
           c("feature", "group", "pct_in", "pct_out"),
           c("symbol", "cluster", "perc", "perc.ctl"))
  stats[, cluster:= factor(cluster, levels(Idents(dat)))]
  stats[, perc.diff:= perc-perc.ctl]

  # Remove unwanted clusters ----
  if(!is.null(select.clusters)) {
    stats <- stats[cluster %in% select.clusters]
    cl <- factor(cl[cl %in% select.clusters], select.clusters)
  }

  # Count active cells ----
  stats$Nact <- unlist(lapply(levels(cl), function(x) {
    Matrix::rowSums(SCT[, cl==x, drop = FALSE]>0)
  }))
  stats$Ninact <- unlist(lapply(levels(cl), function(x) {
    Matrix::rowSums(SCT[, cl==x, drop = FALSE]==0)
  }))
  stats$Nact.ctl <- unlist(lapply(levels(cl), function(x) {
    Matrix::rowSums(SCT[, cl!=x, drop = FALSE]>0)
  }))
  stats$Ninact.ctl <- unlist(lapply(levels(cl), function(x) {
    Matrix::rowSums(SCT[, cl!=x, drop = FALSE]==0)
  }))
  stats$Nclust <- rep(table(cl), each= nrow(SCT))
  stats$Ntotal <- ncol(SCT)

  # Add fisher test ----
  if(add.fisher) {
    stats[, c("OR", "pval"):= {
      fisher.test(
        matrix(
          c(Nact, Ninact, Nact.ctl, Ninact.ctl),
          ncol= 2,
          byrow = T
        )+1
      )[c("estimate", "p.value")]
    }, .(Nact, Ninact, Nact.ctl, Ninact.ctl)]
    stats[, log2OR:= log2(OR)]
    stats[, padj.fisher:= p.adjust(pval, method = "fdr")]
    stats$OR <- stats$pval <- NULL
  }

  # Bind and order ----
  setcolorder(
    stats,
    if(add.fisher)
      c("symbol", "cluster", "avgExpr", "logFC", "padj", "auc", "log2OR", "padj.fisher") else
        c("symbol", "cluster", "avgExpr", "logFC", "padj", "auc")
  )

  # Return
  return(stats)
}
