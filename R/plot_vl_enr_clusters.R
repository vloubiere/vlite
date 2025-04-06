#' @describeIn vl_motif_cl_enrich method to plot enrichment objects (containing variable, log2OR and padj)
#' @export
plot.vl_enr_cl <- function(obj,
                           x.breaks,
                           padj.cutoff= 0.05,
                           top.enrich= Inf,
                           min.counts= 3L,
                           order= "padj",
                           color.breaks,
                           cex.balloons= 1,
                           col= c("blue", "red"),
                           main= NA,
                           plot.empty.clusters= T)
{
  # Checks
  if(!(order %in% c("padj", "log2OR")))
    stop("Possible values for order are 'padj', 'log2OR'")

  # Import and select based on padj and mi counts
  DT <- data.table::copy(obj)
  if(!is.factor(DT$cl))
    DT[, cl:= factor(cl)]
  DT <- DT[set_hit>=min.counts & padj <= padj.cutoff & log2OR > 0]
  # Handle infinite (negative are irrelevant with this plot)
  if(any(is.infinite(DT$log2OR)))
    warning("Attempt to cap infinite log2OR values max/min finite log2OR.
            If plot fails, try another representation or cap manually")
  if(any(DT$log2OR==Inf) && nrow(DT[log2OR>0 & is.finite(log2OR)]))
    DT[log2OR==Inf, log2OR:= max(DT[log2OR>0 & is.finite(log2OR), log2OR])]
  # Check if any enrichment
  if(!nrow(DT))
    stop("No enrichment found with current cutoffs!")
  # Order
  if(order=="padj")
    setorderv(DT, c("cl", "padj")) else if(order=="log2OR")
      DT <- DT[order(cl, -abs(log2OR))]
  # select top.enrich
  if(any(DT[, .N, cl]$N>top.enrich))
    DT <- DT[variable %in% DT[rowid(DT$cl)<=top.enrich, variable]]
  # Save ordering before dcast
  DT[, variable:= factor(variable, levels= unique(variable))]
  # Remove empty clusters
  if(!plot.empty.clusters)
    DT[, cl:= droplevels(cl)]
  # Add y coordinates to DT and return (matrix upside down)
  DT[, y:= max(as.numeric(variable))-as.numeric(variable)+1]
  # dcast
  x <- dcast(DT, variable~cl, value.var = "log2OR", drop= F)
  x <- as.matrix(x, 1)
  rownames(x) <- DT[rownames(x), name, on= "variable", mult= "first"]
  color.var <- dcast(DT, variable~cl, value.var = "padj", drop= F)
  color.var <- as.matrix(color.var, 1)
  color.var <- -log10(color.var)
  # Plot
  balloons_plot(var1 = x,
                var2 = color.var,
                legend.size.breaks = x.breaks,
                col= col,
                cex= cex.balloons,
                main= main,
                var1.legend.title= "OR (log2)",
                var2.legend.title= "padj (-log10)")
  # Return
  invisible(DT)
}
