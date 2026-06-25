#' Volcano Plot
#' 
#' A wrapper around ggrepelScatterplot to plot a nice volcano Plot where specific
#' points can be labelled using ggrepel.
#'
#' @param log2FC 
#' @param padj 
#' @param labels 
#' @param label.col 
#' @param main 
#' @param log2FC.cutoff 
#' @param padj.cutoff 
#'
#' @returns A volcano plot where specific points can be labelled
#' @export
#'
#' @examples
volcanoPlot <- function(
    log2FC,
    padj,
    labels = rep("", length(log2FC)),
    label.col = "black",
    main= NULL,
    log2FC.cutoff= 1,
    padj.cutoff= 0.05
) {
  # Cap log2FC and padj
  if(any(log2FC==Inf | log2FC==(-Inf))) {
    warning("Capping infinite log2FC values!")
    log2FC[log2FC==Inf] <- max(c(log2FC.cutoff, log2FC[is.finite(log2FC)]), na.rm= T)
    log2FC[log2FC==(-Inf)] <- min(c(-log2FC.cutoff, log2FC[is.finite(log2FC)]), na.rm= T)
  }
  log10.padj <- -log10(padj)
  if(any(log10.padj==Inf | log10.padj==(-Inf))) {
    warning("Capping infinite log10.padj values!")
    log10.padj[log10.padj==Inf] <- max(c(padj.cutoff, log10.padj[is.finite(log10.padj)]), na.rm= T)
    log10.padj[log10.padj==(-Inf)] <- min(c(-padj.cutoff, log10.padj[is.finite(log10.padj)]), na.rm= T)
  }
  
  # Compute color
  col <- ifelse(padj <= padj.cutoff & abs(log2FC) >= log2FC.cutoff, "tomato", "lightgrey")
  
  # Plot
  vlite::ggrepelScatterplot(
    log2FC,
    log10.padj,
    col= col,
    xlab= "Fold change (log2)",
    ylab= "FDR (-log10)",
    main= main,
    label = labels,
    label.col = label.col
  )
}