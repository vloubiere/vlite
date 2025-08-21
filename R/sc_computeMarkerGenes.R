#' Compute Marker Genes for Selected Clusters
#'
#' Identifies marker genes for specified clusters in a Seurat object using Wilcoxon rank-sum test (via presto) and, optionally, Fisher's exact test.
#' Performs "cluster vs rest" comparisons.
#'
#' @param dat A Seurat object (must contain an SCT assay).
#' @param select.clusters Character vector of cluster names to analyze. If missing, marker genes are computed for all clusters, by comparing cells
#' within each cluster to the cells outside (except cells from excluded clusters, see next argument).
#' @param exclude.clusters Specified clusters will be excluded before any analysis.
#' @param pairwise Should selected clusters be merged and compared to non-cluster cells using a single pairwise comparison? Default= FALSE.
#' If a character vector is provided, those clusters will be used as the control group.
#' @param value.var A factor vector of clusters of length ncol(dat). Default= Idents(dat).
#' @param add.fisher Logical; if TRUE, adds Fisher's exact test statistics for each gene/cluster.
#' @param verbose For debugging purposes.
#'
#' @return A data.table with marker gene statistics for each cluster, including average expression, log fold change, adjusted p-values, AUC, and
#' optionally Fisher's test results.
#'
#' @details
#' - Uses the SCT assay's data slot for expression values.
#' - Only clusters in `select.clusters` are analyzed if provided.
#' - Fisher's test is performed on active/inactive cell counts if `add.fisher = TRUE`.
#'
#' @example
#' # Import
#' dat <- readRDS("db/seurat_individual_samples/final/PH18_final.rds")
#'
#' # Find top markers for all clusters
#' sc_computeMarkerGenes(dat)
#'
#' # For a subset of cluster
#' sc_computeMarkerGenes(dat, select.clusters = c("a3-r", "PPN"))
#'
#' # Using pairwise comparison
#' sc_computeMarkerGenes(dat, select.clusters = c("a3-r", "PPN"))
#'
#' @export
sc_computeMarkerGenes <- function(
    dat,
    select.clusters= NULL,
    exclude.clusters= NULL,
    pairwise= FALSE,
    value.var= Idents(dat),
    add.fisher= FALSE,
    cleanup.cache= FALSE,
    verbose= FALSE
)
{
  # Checks ----
  stopifnot(is.logical(add.fisher))
  if(class(dat)!="Seurat")
    stop("dat should be a seurat object.")
  if(!is.factor(value.var))
    value.var <- factor(value.var, unique(value.var))
  # Coerce clusters to characters
  if(!is.null(select.clusters))
    select.clusters <- as.character(select.clusters)
  if(!is.null(exclude.clusters))
    exclude.clusters <- as.character(exclude.clusters)
  if(isTRUE(pairwise))
    pairwise <- setdiff(levels(value.var), c(select.clusters, exclude.clusters))
  if(!isFALSE(pairwise))
    pairwise <- as.character(pairwise)
  if(isFALSE(pairwise))
    pairwise <- NULL
  # Check they exist & do not overlap
  if(!all(c(select.clusters, exclude.clusters, pairwise) %in% levels(value.var)))
    stop("select.clusters, exclude.clusters and pairwise clusters should all match value.var level(s).")
  if(anyDuplicated(c(select.clusters, exclude.clusters, pairwise)))
    stop("Unwanted overlap(s) between select.clusters, exclude.clusters and pairwise clusters.")

  # Use a temp directory for caching ----
  set.params <- list(digest::digest(dat@assays$SCT$data),
                     select.clusters,
                     exclude.clusters,
                     pairwise,
                     value.var,
                     add.fisher)
  set.key <- digest::digest(set.params)
  set.cache <- file.path(tempdir(), paste0(set.key, ".rds"))

  # Compute marker genes ----
  if(cleanup.cache | !file.exists(set.cache)) {

    # Import data ----
    SCT <- dat@assays$SCT$data
    if(nrow(SCT)==1)
      stop("SCT contains only one row.")

    # Remove excluded clusters ----
    if(!is.null(exclude.clusters)) {
      SCT <- SCT[, !value.var %in% exclude.clusters]
      if(verbose)
        print(table(value.var))
      value.var <- value.var[!value.var %in% exclude.clusters]
      value.var <- droplevels(value.var)
      if(verbose)
        print(table(value.var))
    }

    # Simplify clusters and data based on selection ----
    if(!is.null(pairwise)) {
      # Subset SCT and clusters
      SCT <- SCT[, value.var %in% c(select.clusters, pairwise)]
      value.var <- value.var[value.var %in% c(select.clusters, pairwise)]
      # Unique selection name
      selection.name <- ifelse(length(select.clusters)==1, select.clusters, "selection")
      # Simplify clusters
      value.var <- ifelse(value.var %in% select.clusters, selection.name, "control")
      # Select.cluster can now be set to selection.name
      select.clusters <- selection.name
      # As factor
      value.var <- factor(value.var, c(selection.name, "control"))
    } else if(!is.null(select.clusters)) {
      # Unique control name
      ctl.name <- rev(make.unique(c(select.clusters, "control")))[1]
      # Simplify clusters
      value.var <- ifelse(value.var %in% select.clusters, as.character(value.var), ctl.name)
      # As factor
      value.var <- factor(value.var, c(select.clusters, ctl.name))
    }
    if(verbose)
      print(table(value.var))

    # Compute avgExpr, wilcox and auc ----
    stats <- presto::wilcoxauc(SCT, y= value.var)
    stats <- as.data.table(stats)
    stats$statistic <- stats$pval <- NULL
    setnames(stats,
             c("feature", "group", "pct_in", "pct_out"),
             c("symbol", "cluster", "perc", "perc.ctl"))
    stats[, cluster:= factor(cluster, levels(value.var))]
    stats[, perc.diff:= perc-perc.ctl]

    # Remove unwanted clusters ----
    if(!is.null(select.clusters))
      stats <- stats[cluster %in% select.clusters]

    # Count active cells ----
    stats[, c("Nact", "Ninact", "Nact.ctl", "Ninact.ctl"):= {
      idx <- as.character(value.var)==as.character(cluster)
      .(Matrix::rowSums(SCT[, idx, drop = FALSE]>0),
        Matrix::rowSums(SCT[, idx, drop = FALSE]==0),
        Matrix::rowSums(SCT[, !idx, drop = FALSE]>0),
        Matrix::rowSums(SCT[, !idx, drop = FALSE]==0))
    }, cluster]
    stats[, Nclust:= Nact+Ninact]

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

    # Save
    saveRDS(stats, set.cache)
  } else
    stats <- readRDS(set.cache)

  # Return
  invisible(stats)
}
