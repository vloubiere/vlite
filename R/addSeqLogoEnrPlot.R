#' Plot motif logos on the left of an enrichment plot
#'
#' Specialized funciton used to plot motif logs on the left of a heatmap/balloon plot.
#'
#' @param obj Either a data.table containing 'name' and 'y' columns, similar to the $rows data.table
#' returned by vl_heatmap or the `vl_enr_cl` object returned by the ?plot.vl_enr_cl function.
#' The 'name' values should match the motif names present in the provided ICMatrixList (see next arugment),
#' and y corresponds to the y plotting positions of each row.
#' @param ICMatrixList An ICMatrixList object containing the motifs to add.
#' @param cex.width A width expansion factor applied to letter heights.
#' @param cex.height A height expansion factor applied to letter heights.
#'
#' @examples
#' # Select br PPM
#' sel <- plot(
#'   enr,
#'   top.enrich = 5,
#'   order = "log2OR",
#'   padj.cutoff = 1e-3,
#'   col= c("pink", "red"),
#'   cex = .4,
#'   main = "ChIP motif enrich."
#' )
#' addSeqLogoEnrPlot(
#'   obj = sel,
#'   ICM = mot$ICM
#' )
#' 
#' @export
addSeqLogoEnrPlot <- function(
    obj,
    ICM,
    cex.width= 1,
    cex.height= 1
)
{
  # Checks ----
  stopifnot(is.data.table(obj))
  obj <- data.table::copy(obj)
  if("y.pos" %in% names(obj) && !"y" %in% names(obj))
    setnames(obj, "y.pos", "y")
  stopifnot(c("name", "y") %in% names(obj))
  stopifnot(inherits(ICM, "ICMatrixList"))
  
  # Plot ----
  obj[, {
    addSeqLogo(
      mat = ICM[[name]],
      x = par("usr")[1]-strwidth(paste0(name, "M"), cex= par("cex.axis"))-diff(grconvertX(c(0, par("mgp")[2]), "line", "user")),
      y = y,
      pos= 2,
      cex.width = .6*cex.width,
      cex.height = .6*cex.height
    )
  }, .(name, y)]
}
