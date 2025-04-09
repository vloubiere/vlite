#' Plot Clustered Enrichment Results for vl_enr_cl Objects
#'
#' This function generates a balloon plot to visualize clustered enrichment results from an object of class `vl_enr_cl`.
#' The balloons represent the log2 odds ratio (log2OR) of enrichment, with their size corresponding to the odds ratio
#' and their color representing the adjusted p-values (-log10 transformed).
#'
#' @param obj An object of class `vl_enr_cl`. This should be a data table containing clustered enrichment results,
#' including columns for `padj` (adjusted p-values), `log2OR` (log2 odds ratio), `set_hit` (counts),
#' `name` (feature names), and `cl` (clusters).
#' @param padj.cutoff Numeric. The adjusted p-value cutoff to filter enrichments before plotting.
#' Default is `0.05`.
#' @param top.enrich Integer. The maximum number of top enrichments to plot per cluster, based on the `order` parameter.
#' Default is `Inf` (plot all enrichments that pass the cutoffs).
#' @param min.counts Integer. The minimum number of counts (`set_hit`) required to include an enrichment
#' in the plot. Default is `3L`.
#' @param order Character. The metric used to order enrichments before selecting the top enrichments.
#' Possible values are `padj` (adjusted p-value) or `log2OR` (log2 odds ratio). Default is `log2OR`.
#' @param color.breaks Numeric vector. The color breaks to use for the color scale. Defaults to the range of
#' filtered `padj` values (-log10 transformed).
#' @param size.legend.breaks Numeric vector. Breakpoints for the size legend in the balloon plot. Default is `NULL`.
#' @param cex Numeric. Scaling factor for the plot. Default is `1`.
#' @param col Character vector. The color scale to use for the balloon plot. Default is `c("blue", "red")`.
#' @param main Character. The main title for the plot. Default is `NA` (no title).
#' @param plot.empty.clusters Logical. Whether to include clusters with no enrichments in the plot.
#' Default is `TRUE`.
#'
#' @details
#' The function filters the input data based on the specified `padj.cutoff`, `min.counts`, and `log2OR > 0` thresholds.
#' It then orders the data based on the `order` parameter and selects the top enrichments per cluster.
#' The balloon plot displays the log2 odds ratio (log2OR) as the size of the balloons, with colors representing
#' the adjusted p-values (-log10 transformed). Clusters with no enrichments can be optionally excluded from the plot.
#'
#' @return Invisibly returns the filtered and plotted data table.
#'
#' @examples
#' # Example usage:
#' For an example using GO enrichment, see ?vl_GOenrich()
#'
#' # Assuming `clustered_enrichment_results` is an object of class `vl_enr_cl`:
#' plot.vl_enr_cl(clustered_enrichment_results, padj.cutoff = 0.01, top.enrich = 5, min.counts = 5)
#'
#' @export
plot.vl_enr_cl <- function(obj,
                           padj.cutoff= 0.05,
                           top.enrich= Inf,
                           min.counts= 3L,
                           order= "log2OR",
                           color.breaks= NULL,
                           size.legend.breaks= NULL,
                           cex= 1,
                           col= c("blue", "red"),
                           main= NA,
                           plot.empty.clusters= T)
{
  # Checks
  if(!(order %in% c("padj", "log2OR")))
    stop("Possible values for order are 'padj', 'log2OR'")
  if(uniqueN(obj$name)!=uniqueN(obj$variable))
    stop("Several 'variable' values have similar 'name' value in obj. 'name' values should be unique!")

  # Import and select based on padj and min.counts
  DT <- data.table::copy(obj)
  DT <- DT[set_hit>=min.counts & padj <= padj.cutoff & log2OR > 0]

  # Checks
  if(any(is.infinite(DT$log2OR)))
    stop("Infinite enrichment values should be capped before plotting.")
  if(!nrow(DT))
    stop("No enrichment found with current cutoffs!")

  # Order
  if(order=="padj"){
    setorderv(DT, c("cl", "padj"))
  } else if(order=="log2OR") {
    DT <- DT[order(cl, -abs(log2OR))]
  }

  # select top.enrich
  if(any(DT[, .N, cl]$N>top.enrich))
    DT <- DT[name %in% DT[rowid(DT$cl)<=top.enrich, name]]

  # Save ordering before dcast
  DT[, name:= factor(name, levels= unique(name))]

  # Remove empty clusters
  if(!plot.empty.clusters)
    DT[, cl:= droplevels(cl)]

  # Add y coordinates to DT and return (matrix upside down)
  DT[, y:= max(as.numeric(name))-as.numeric(name)+1]

  # Dcast
  log2OR <- dcast(DT, name~cl, value.var = "log2OR", drop= F)
  log2OR <- as.matrix(log2OR, 1)
  padj <- dcast(DT, name~cl, value.var = "padj", drop= F)
  padj <- as.matrix(padj, 1)
  padj <- -log10(padj)

  # Plot
  balloons_plot(size.var = log2OR,
                color.var = padj,
                color.breaks = color.breaks,
                size.legend.title = "OR (log2)",
                color.legend.title = "padj (-log10)",
                size.legend.breaks = size.legend.breaks,
                cex= cex,
                main= main)

  # Return
  invisible(DT)
}
