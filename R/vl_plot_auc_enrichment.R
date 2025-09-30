#' Title
#'
#' @param auc A vl_auc object created by ?vl_motifAUC.
#' @param cutoff Cutoff on the maximum value across all conditions.
#' @param min.diff.cutoff Cutoff on the minimum difference between lowest and highest value.
#' @param topN After sorting the table using order.var (decreasing= T), subset only topN values for ploting.
#' @param order.var Value used to sort the data before subsetting topN values. Possible values are 'score' (use
#' the absolute of the value.var) or 'diff' (sort by difference between highest and lowest value.var).
#' @param value.var The value to be plotted. Possible value are "NES" or "auc". Default= "NES".
#' @param ... Extra arguments to be passed to the heatmap funciton.
#'
#' @return Invisible returns the final matrix.
#' @export
#'
#' @examples
plot.vl_auc <- function(auc,
                        cutoff= ifelse(value.var=="NES", 3, .6),
                        min.diff.cutoff= 0,
                        topN= Inf,
                        order.var= "score",
                        value.var= "NES",
                        cluster.rows= T,
                        cluster.cols= F,
                        ranks.as.rows= FALSE,
                        kmeans.k= NA,
                        show.row.clusters= "right",
                        show.numbers= NULL,
                        legend.cex = 1,
                        numbers.cex = .5,
                        pdf.file= NULL)
{
  # Checks ----
  stopifnot(value.var %in% c("auc", "NES"))
  stopifnot(order.var %in% c("score", "diff"))

  # Format ----
  cols <- c(names(auc)[1:2], value.var)
  auc <- auc[, cols, with= F]
  if(ranks.as.rows)
    setnames(auc, c("var", "rank", "score")) else
      setnames(auc, c("rank", "var", "score"))

  # Apply cutoffs ----
  auc[, diff:= diff(range(score)), var]
  sel1 <- auc[abs(score) >= cutoff, var]
  sel2 <- auc[, diff >= min.diff.cutoff, var][(V1), var]
  auc <- auc[var %in% intersect(sel1, sel2)]

  # Select top values ----
  if(topN != Inf) {
    auc <- if(order.var == "score")
      auc[order(abs(score), decreasing= T)] else
        auc[order(diff, decreasing= T)]
    top <- auc[, .SD[1:min(c(.N, topN))], rank]$var
    auc <- auc[var %in% top]
  }

  # Select lines base on minimum value and range ----
  if(!nrow(auc))
    stop("No enrichment found with actual cutoff.")

  # Dcast ----
  mat <- dcast(auc, var~rank, value.var = "score")
  mat <- as.matrix(mat, 1)
  print(paste(nrow(mat), " enrichments found."))

  # Plotting variable ----
  breaks <- max(abs(mat), na.rm= T)
  breaks <- ifelse(breaks>15, 15, ifelse(breaks<5, 5, breaks))
  breaks <- c(-breaks, -3, 3, breaks)
  Cc <- circlize::colorRamp2(breaks, c("blue", "white", "white", "red"))
  breaks <- seq(breaks[1], rev(breaks)[1], length.out= 21)
  col <- Cc(breaks)

  # Add numbers ----
  if(is.null(show.numbers))
    show.numbers <- round(mat)

  # Plot ----
  hm <- vl_heatmap(mat,
                   legend.title = value.var,
                   legend.cex = legend.cex,
                   show.numbers = show.numbers,
                   numbers.cex= numbers.cex,
                   breaks= breaks,
                   col= col,
                   kmeans.k = kmeans.k,
                   cluster.rows = cluster.rows,
                   cluster.cols = cluster.cols,
                   show.legend = "top",
                   pdf.file = pdf.file,
                   show.row.clusters = show.row.clusters)

  # Return
  invisible(mat)
}
