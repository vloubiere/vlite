#' Subset top marker genes
#'
#' Subset top marker genes from a data.table generate by >sc_computeMarkerGenes.
#'
#' @param dat A data.table with marker gene statistics. See ?sc_computeMarkerGenes.
#' @param topN Number of top genes per cluster. Default= Inf (all markers).
#' @param perc.cutoff Minimum percentage of expressing within cluster.
#' @param perc.diff.cutoff Minimum difference in percentage of expressing cells within cluster vs. outside cluster.
#' @param padj.wilcox.cutoff Adjusted Wilcoxon p-value cutoff. Default= 0.05.
#' @param logFC.cutoff Log fold change cutoff. Default= 0.25.
#' @param padj.fisher.cutoff Adjusted Fisher p-value cutoff. Default= NULL.
#' @param log2OR.cutoff Log2 odds ratio cutoff. Default= NULL.
#' @param auc.cutoff AUC cutoff. Default= NULL.
#' @param order.var Variable to order genes by ("perc.diff", "auc", "logFC", "log2OR").
#' @param alternative The direction in which filtering should be done. SHould be one of 'greater' or 'smaller'.
#' @return The subset of the rds object corresponding to top marker genes.
#' @export

sc_topMarkers <- function(
    dat,
    topN= Inf,
    perc.cutoff= 10,
    perc.diff.cutoff= 20,
    padj.wilcox.cutoff= 0.05,
    logFC.cutoff= .25,
    padj.fisher.cutoff= NULL,
    log2OR.cutoff= NULL,
    auc.cutoff= NULL,
    order.var= "perc.diff",
    alternative= "greater"
) {
  # Checks ----
  if(!is.data.table(dat))
    stop("dat should be a data.table containing marker genes.")
  missing <- setdiff(
    c(order.var, "padj",  "logFC", "perc", "perc.diff", "auc"),
    names(dat)
  )
  if(length(missing))
    stop(paste0("Missing columns: ", paste0(missing, collapse= ", ")))
  if(!order.var %in% c("perc.diff", "auc", "logFC", "log2OR"))
    stop("order.var should be one of: perc.diff, auc, logFC, log2OR")
  if(!alternative %in% c("greater", "smaller"))
    stop("Alternative should be one of 'greater' or 'smaller'")
  if(!is.null(logFC.cutoff) && logFC.cutoff<0)
    stop("logFC.cutoff should be >=0 (automatically set to negative if alternative= 'smaller')")
  if(!is.null(log2OR.cutoff) && log2OR.cutoff<0)
    stop("log2OR.cutoff should be >=0 (automatically set to negative if alternative= 'smaller')")
  if(!is.null(perc.cutoff) && !between(perc.cutoff, 0, 100))
    stop("perc.cutoff should be a numeric value between 0 and 100 (automatically set to 100-perc.cutoff if alternative= 'smaller')")
  if(!is.null(perc.diff.cutoff) && !between(perc.diff.cutoff, 0, 100))
    stop("perc.diff.cutoff should be a numeric value between 0 and 100 (automatically set to negative if alternative= 'smaller')")
  if(!is.null(padj.wilcox.cutoff) && !between(padj.wilcox.cutoff, 0, 1))
    stop("padj.wilcox.cutoff should be a numeric value between 0 and 1.")
  if(!is.null(padj.fisher.cutoff) && !between(padj.fisher.cutoff, 0, 1))
    stop("padj.fisher.cutoff should be a numeric value between 0 and 1.")
  if((any(!is.null(c(log2OR.cutoff, padj.fisher.cutoff))) && !all(c("padj.fisher", "log2OR") %in% names(dat))))
    stop("padj.fisher/log2OR column(s) missing -> cutoff cannot be applied.")
  if(!is.null(auc.cutoff) && !between(auc.cutoff, 0, 1))
    stop("auc.cutoff should be a numeric value between 0 and 1 (automatically set to 1-auc.cutoff if alternative= 'smaller').")

  # Import data ----
  dat <- data.table::copy(dat)

  # Order ----
  setorderv(dat, order.var, ifelse(alternative=="smaller", 1, -1))
  dat[, rank:= seq(.N), cluster]

  # Copy ----
  top <- data.table::copy(dat)

  # Successive cutoffs ----
  top <- if(topN==0)
    top[0] else {

    # Cutoff on wilcox test (expression) ----
    if(!is.null(padj.wilcox.cutoff))
      top <- top[padj <= padj.wilcox.cutoff]
    if(!is.null(logFC.cutoff)) {
      top <- if(alternative=="greater")
        top[logFC >= logFC.cutoff] else
          top[logFC <= (-logFC.cutoff)]
    }

    # Cutoff on fisher test (over-representation expressing cells) ----
    if(!is.null(padj.fisher.cutoff))
      top <- top[padj.fisher <= padj.fisher.cutoff]
    if(!is.null(log2OR.cutoff)) {
      top <- if(alternative=="greater")
        top[log2OR >= log2OR.cutoff] else
          top[log2OR <= (-log2OR.cutoff)]
    }

    # Cutoff on percentage of expressing cells ----
    if(!is.null(perc.cutoff)) {
      top <- if(alternative=="greater")
        top[perc >= perc.cutoff] else
          top[perc <= (100-perc.cutoff)]
    }

    # Cutoff on difference in percentage of expressing cells ----
    if(!is.null(perc.diff.cutoff)) {
      top <- if(alternative=="greater")
        top[perc.diff >= perc.diff.cutoff] else
          top[perc.diff <= (-perc.diff.cutoff)]
    }

    # AUC cutoff ----
    if(!is.null(auc.cutoff)) {
      top <- if(alternative=="greater")
        top[auc >= auc.cutoff] else
          top[auc <= 1-auc.cutoff]
    }

    # Select N top genes ----
    if(nrow(top))
      top <- top[, .SD[1:min(c(.N, topN))], keyby= cluster]
  }

  # Return ----
  return(top)
}
