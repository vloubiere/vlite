#' Plot Clustered Enrichment Results for vl_enr_cl Objects
#'
#' This function generates a balloon plot to visualize clustered enrichment results from an object of class 'vl_enr_cl'.
#' The balloons represent the log2 odds ratio (log2OR) of enrichment, with their size corresponding to the odds ratio
#' and their color representing the adjusted p-values (-log10 transformed).
#'
#' @param obj A data.table of class 'vl_enr_cl' containing columns 'padj', 'log2OR', 'set_hit', 'name' and 'cl'.
#' @param log2OR.cutoff Numeric cutoff used to filter positive log2OR values before plotting.
#' @param padj.cutoff Numeric cutoff used to filter adjusted p-values before plotting.
#' Default= 0.05.
#' @param top.enrich Integer specifying the maximum number of top enrichments to plot, based on the
#' 'order' parameter. Default= Inf (no selection).
#' @param min.counts Integer speifying the minimum number of counts required to include an enrichment
#' in the plot ('set_hit' column). Default= 3.
#' @param order The metric that should be used to order enrichments before selecting the top enrichments.
#' Possible values are 'padj' or 'log2OR'. Default= 'log2OR'.
#' @param color.breaks Color breaks to use for the color scale. Defaults to the range of filtered
#' -log10(adjusted p-values).
#' @param col The color scale to use for the balloon plot. Default= c('blue', 'red').
#' @param size.legend.breaks Breakpoints for the size legend in the balloon plot. Defaults to the range
#' of filtered log2OR values.
#' @param cex Numeric expansion factor used to adjust balloons' sizes. Default= 1.
#' @param main Character. The main title for the plot. Default= NA.
#' @param plot.empty.clusters Should empty clusters be plotted? Default= TRUE.
#'
#' @return Invisibly returns the filtered and plotted data table.
#'
#' @examples
#' # Example usage:
#' For an example using GO enrichment, see ?vl_GOenrich()
#' For an example using motif enrichment, see ?vl_motifEnrich()
#'
#' # Assuming `clustered_enrichment_results` is an object of class `vl_enr_cl`:
#' plot.vl_enr_cl(clustered_enrichment_results, padj.cutoff = 0.01, top.enrich = 5, min.counts = 5)
#'
#' @export
plot.vl_enr_cl <- function(obj,
                           log2OR.cutoff= 0,
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
  if(any(obj[, .N, .(name, cl)]$N > 1))
    stop("Several lines were found with similar name / cl combination.")
  if(any(is.na(obj$name)))
    stop("Some names in DT are NA.")

  # Import and select based on padj and min.counts
  DT <- data.table::copy(obj)
  DT <- DT[set_hit >= min.counts & padj <= padj.cutoff & log2OR >= log2OR.cutoff]

  # Checks
  if(any(is.infinite(DT$log2OR)))
    stop("Infinite enrichment values should be capped before plotting.")
  if(any(DT$padj==0))
    stop("Some padjust are equal to 0 and should be set to a minimum positive value before plotting.")
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
                col = col,
                cex= cex,
                main= main)

  # Return
  invisible(DT)
}
