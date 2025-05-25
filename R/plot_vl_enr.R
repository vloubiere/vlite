#' Plot Enrichment Results for vl_enr Objects
#'
#' This function generates a barplot to visualize enrichment results from an object of class 'vl_enr'.
#' The bars represent the log2 odds ratio (log2OR) of enrichment, and the color scale corresponds to
#' the adjusted p-values (-log10 transformed).
#'
#' @param obj A data.table of class 'vl_enr' containing columns 'padj', 'log2OR' and 'set_hit'.
#' See ?vl_motifEnrich for an example.
#' @param log2OR.abs.cutoff Numeric cutoff used to filter absolute log2OR values before plotting.
#' @param padj.cutoff Numeric cutoff used to filter adjusted p-values before plotting.
#' Default= 0.05.
#' @param top.enrich Integer specifying the maximum number of top enrichments to plot, based on the
#' 'order' parameter. Default= Inf (no selection).
#' @param min.counts Integer speifying the minimum number of counts required to include an enrichment
#' in the plot ('set_hit' column). Default= 3.
#' @param order The metric that should be used to order enrichments before selecting the top enrichments.
#' Possible values are 'padj' or 'log2OR'. Default= 'log2OR'.
#' @param xlab Label for the x-axis of the barplot. Default= 'Odd Ratio (log2)'.
#' @param col Color scale to use for the barplot. Default= c('blue', 'red').
#' @param breaks Numeric vector specifying the color breaks to use for the color scale. Defaults to the range of
#' filtered -log10(adjusted p-values).
#'
#' @return Invisibly returns the filtered and plotted data table.
#'
#' @examples
#' # Example usage:
#' For an example using GO enrichment, see ?vl_GOenrich()
#' For an example using motif enrichment, see ?vl_motifEnrich()
#'
#' # Assuming `enrichment_results` is an object of class `vl_enr`:
#' plot.vl_enr(enrichment_results, padj.cutoff = 0.01, top.enrich = 10, min.counts = 5)
#'
#' @export
plot.vl_enr <- function(obj,
                        log2OR.abs.cutoff= 0,
                        padj.cutoff= 0.05,
                        top.enrich= Inf,
                        min.counts= 3L,
                        order= "log2OR",
                        xlab= "Odd Ratio (log2)",
                        breaks= NULL,
                        col= c("blue", "red"))
{
  # Checks
  if(!(order %in% c("padj", "log2OR")))
    stop("Possible values for order are 'padj', 'log2OR'")
  if(any(obj[, .N, .(name, cl)]$N > 1))
    stop("Several lines were found with similar name.")

  # Import and select based on padj and min.counts cutoff
  DT <- data.table::copy(obj)
  DT <- DT[padj<=padj.cutoff & set_hit>=min.counts & abs(log2OR) >= abs(log2OR.abs.cutoff)]

  # Checks
  if(any(is.infinite(DT$log2OR)))
    stop("Infinite enrichment values should be capped before plotting.")
  if(!nrow(DT))
    stop("No enrichment found with current cutoffs!")

  # Order
  if(order=="padj") {
    setorderv(DT, "padj")
  } else if(order=="log2OR") {
    DT <- DT[order(-abs(log2OR))]
  }

  # Select top.enrich
  if(nrow(DT)>top.enrich)
    DT <- DT[seq(nrow(DT))<=top.enrich]

  # Plot
  if(is.null(breaks))
  {
    breaks <- range(-log10(DT$padj), na.rm= T)
    if(length(unique(breaks))==1)
      breaks <- breaks+c(-0.5,0.5)
    breaks <- seq(min(breaks),
                  max(breaks),
                  length.out= length(col))
  }
  Cc <- circlize::colorRamp2(breaks, colorRampPalette(col)(length(breaks)))

  # Reorder by Log2OR before plotting
  setorderv(DT, "log2OR")

  # Barplot
  DT[, y:= barplot(log2OR,
                   horiz= T,
                   names.arg= name,
                   border= NA,
                   col= Cc(-log10(padj)),
                   las= 1,
                   xlab= xlab)]

  # Plot heatkey
  hk.breaks <- seq(breaks[1], breaks[2], length.out= 100)
  heatkey(breaks = hk.breaks,
          col = Cc(hk.breaks))

  # Return
  invisible(DT)
}
