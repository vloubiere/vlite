#' Compute Pariwise Marker Genes
#'
#' Identifies marker genes for a merged group of clusters.
#'
#' @param dat A Seurat object (must contain an SCT assay).
#' @param value.var A factor vector of clusters of length ncol(dat). Default= Idents(dat).
#' @param select.clusters Character vector of cluster names to analyze. If missing, marker genes are computed for all clusters, by
#' comparing observations within all the selected clusters vs. outside.
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
sc_computePariwiseMarkerGenes <- function(
    dat,
    select.clusters,
    ctl.clusters= NULL,
    value.var= Idents(dat),
    add.fisher= FALSE
)
{
  # Checks
  if(class(dat)!="Seurat")
    stop("dat should be a seurat object.")
  if(!is.factor(value.var))
    value.var <- factor(value.var, unique(value.var))
  select.clusters <- unique(select.clusters)
  if(!is.character(select.clusters) || !all(select.clusters %in% levels(value.var)))
    stop("select.clusters should be a character vector matching any level in value.var.")
  if(is.null(ctl.clusters))
    ctl.clusters <- setdiff(levels(value.var), select.clusters)
  if(!length(ctl.clusters))
    stop("At least one control cluster should be used as control")
  if(any(select.clusters %in% ctl.clusters))
    stop("A cluster level was specified both in the selected and control set.")
  print(paste("Tested clusters:", paste0(select.clusters, collapse = ", ")))
  print(paste("Control clusters:", paste0(ctl.clusters, collapse = ", ")))

  # Import data ----
  SCT <- dat@assays$SCT$data
  if(nrow(SCT)==1)
    stop("SCT contains only one row.")

  # Subset clusters and simplify ----
  cl <- fcase(value.var %in% select.clusters, "selected",
              value.var %in% ctl.clusters, "control",
              default = as.character(NA))
  browser()


  # Simplify when using selected clusters ----
  if(!is.null(select.clusters)) {
    non.selected <- rev(make.unique(c(select.clusters, "non.selected")))[1]
    cl <- ifelse(as.character(value.var) %in% select.clusters, as.character(value.var), non.selected)
    cl <- factor(cl, c(select.clusters, non.selected))
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
