#' Plot Enrichment Results for vl_enr Objects
#'
#' This function generates a barplot to visualize enrichment results from an object of class `vl_enr`.
#' The bars represent the log2 odds ratio (log2OR) of enrichment, and the color scale corresponds to
#' the adjusted p-values (-log10 transformed).
#'
#' @param obj An object of class `vl_enr`. This should be a data table containing enrichment results,
#' including columns for `padj` (adjusted p-values), `log2OR` (log2 odds ratio), and `set_hit` (counts).
#' @param padj.cutoff Numeric. The adjusted p-value cutoff to filter enrichments before plotting.
#' Default is `0.05`.
#' @param top.enrich Integer. The maximum number of top enrichments to plot, based on the `order` parameter.
#' Default is `Inf` (plot all enrichments that pass the cutoffs).
#' @param min.counts Integer. The minimum number of counts (`set_hit`) required to include an enrichment
#' in the plot. Default is `3L`.
#' @param order Character. The metric used to order enrichments before selecting the top enrichments.
#' Possible values are `padj` (adjusted p-value) or `log2OR` (log2 odds ratio). Default is `log2OR`.
#' @param xlab Character. The label for the x-axis of the barplot. Default is `Odd Ratio (log2)`.
#' @param col Character vector. The color scale to use for the barplot. Default is `c("blue", "red")`.
#' @param breaks Numeric vector. The color breaks to use for the color scale. Defaults to the range of
#' filtered `padj` values (-log10 transformed).
#'
#' @details
#' The function filters the input data based on the specified `padj.cutoff` and `min.counts` thresholds.
#' It then orders the data based on the `order` parameter and selects the top enrichments to plot.
#' The barplot displays the log2 odds ratio (log2OR) for each enrichment, with colors representing
#' the adjusted p-values (-log10 transformed). A heat key is also displayed to indicate the color scale.
#'
#' @return Invisibly returns the filtered and plotted data table.
#'
#' @examples
#' # Example usage:
#' For an example using GO enrichment, see ?vl_GOenrich()
#'
#' # Assuming `enrichment_results` is an object of class `vl_enr`:
#' plot.vl_enr(enrichment_results, padj.cutoff = 0.01, top.enrich = 10, min.counts = 5)
#'
#' @export
plot.vl_enr <- function(obj,
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

  # Import and select based on padj and min.counts cutoff
  DT <- data.table::copy(obj)
  DT <- DT[padj<=padj.cutoff & set_hit>=min.counts]

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
          top = DT[.N, y],
          left = par("usr")[2]+strwidth("M"),
          col = Cc(hk.breaks),
          main = "FDR (-log10)")

  # Return
  invisible(DT)
}
