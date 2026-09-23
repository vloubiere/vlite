#' Function method to plot motifEnrichRanks recovery curves
#'
#' @param obj An object as returned by ?motifEnrichRanks.
#' @param sel.motif The motif for which empirical cumulative recovery curve should be plotted.
#' @param col The color for the curve. Defaults= "blue.
#' @param mean.curve.col The color for the curve. Defaults= "grey30".
#' @param sd.curve.col The color for the curve. Defaults= "limegreen".
#'
#' @returns A recovery plot
#' @export
#'
#' @examples
plot.vl_enrRank <- function(
    obj,
    sel.motif,
    sel.group= NULL,
    col= "blue",
    mean.curve.col= "grey30",
    sd.curve.col= "limegreen",
    main= sel.motif
) {

  # Extract vars ----
  R <- obj$recov$motif.recov[,sel.motif]
  mR <- obj$recov$mean.recov
  sR <- obj$recov$mean.recov+(2*obj$recov$sd.recov)
  
  # Compute ylim and initiate plot ----
  ylim <- range(
    c(0, mR, sR, R)
  )
  plot(
    x = c(0, length(R)),
    y= ylim,
    type= "n",
    main= main,
    xlab= "Rank",
    ylab= "Cumulative recovery"
  )
  
  # Plot mean curves ----
  lines(seq(R), y= R, col= col)
  lines(x= seq(mR), y= mR, col= mean.curve.col)
  lines(x= seq(sR), y= sR, col= sd.curve.col)
  
  # Add NES score, FDR and number of put binding sites ----
  NES <- round(obj$enr[motif==sel.motif, NES], 1)
  FDR <- if("padj" %in% names(obj$enr)) 
    formatC(obj$enr[motif==sel.motif, padj], format= "e", digits = 1) else
      NA
  pBS <- obj$enr[motif==sel.motif, leading.edge]
  # plot
  vl_legend(
    "topleft",
    legend= c(
      paste0("NES= ", NES),
      paste0("FDR= ", FDR),
      paste0("pBS= ", pBS)
    ),
    text.col= col
  )
  
  # Add leading edge ----
  leading.edge <- obj$enr[motif==sel.motif, leading.edge]
  if(!is.na(leading.edge)) {
    lines(
      c(0, rep(leading.edge, 2)),
      c(rep(R[leading.edge], 2), 0),
      lty= "33"
    )
  }
}