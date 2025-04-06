#' Compute motif enrichment for several groups
#'
#' Compute motif enrichment for each element of a list (clusters, groups....)
#'
#' @param counts.list List of data.table containing counts for the regions of interest and potentially control regions (see next argument).
#' @param control.cl IDs of clusters to be used as background. Default= NULL, meaning all clusters are used except the one being tested.
#' @param names Convenient names to be used for aggregating and plotting. By default, returns motif_cluster matching motif_ID in vl_Dmel_motifs_DB_full.
#' @param plot Should the result be plotted using balloons plot? Default to FALSE
#' @param padj.cutoff Cutoff for balloons to be plotted.
#' @param top.enrich Select top enriched motifs/cluster. Default= Inf (All).
#' @param min.counts The minimum number of counts required to call a hit. Default= 3L.
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj" (the default) and "log2OR".
#' @param x.breaks Breaks used for balloon's sizes.
#' @param color.breaks Color breaks used for coloring.
#' @param col Vector of colors used for coloring.
#' @param main Title. Default= NA.
#' @param add.motifs Should motif pwms be plotted?
#' @param cex.width Expansion factor for motif widths.
#' @param cex.balloons Expansion factor for balloons.
#' @param plot.empty.clusters Should empty clusters be plotted? Default= TRUE.
#' @param cex.height Expansion factor for motif heights.
#'
#' @examples
#' # Resize example peaks
#' SUHW <- resizeBed(vl_SUHW_top_peaks, genome = "dm3")
#' STARR <- resizeBed(vl_STARR_DSCP_top_peaks, genome = "dm3")
#'
#' # Generate same number of random regions
#' random <- controlRegionsBSgenome(bed= STARR, genome= "dm3")
#'
#' # Combine the three sets of regions
#' combined <- rbindlist(list(SUHW= SUHW,
#'                            STARR= STARR,
#'                            random= random),
#'                       idcol = "cluster",
#'                       fill = T)
#'
#' # Count JAPSPAR motifs
#' counts <- vl_motifCounts(combined, genome= "dm3", sel= pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#'
#' # Motifs can also be counted using a custom PWMatrixList, for example for promoter motifs:
#' prom_db <- readRDS("/groups/stark/almeida/data/motifs/CP_motifs/CP_motifs_PWM.rds")
#' prom <- prom_db$Pwms_log_odds
#' for(i in seq(prom)) # This set needs a fix
#'   prom[[i]]@profileMatrix <- pwmPercToLog(prom_db$Pwms_perc[[i]]@profileMatrix)
#'
#' prom_motifs <- vl_motifCounts(combined,
#'                                pwm_log_odds= prom,
#'                                genome= "dm3")
#'
#' # Enrichment per cluster
#' enr <- vl_motifEnrichClusters(split(counts, combined$cluster),
#'                           control.cl = "random")
#'
#' # Plot
#' vl_par()
#' plot(enr,
#'      padj.cutoff= 1e-5)
#'
#' @return Fisher test data.table.
#' @export
vl_motifEnrichClusters <- function(counts.list,
                                   control.cl= NULL,
                                   names= vl_Dmel_motifs_DB_full[colnames(counts), motif_cluster, on= "motif_ID"],
                                   plot= F,
                                   padj.cutoff= 1e-5,
                                   top.enrich= Inf,
                                   min.counts= 3L,
                                   order= "padj",
                                   x.breaks,
                                   color.breaks,
                                   cex.balloons= 1,
                                   col= c("cornflowerblue", "lightgrey", "tomato"),
                                   main= NA,
                                   plot.empty.clusters= T,
                                   add.motifs= F,
                                   cex.width= 1,
                                   cex.height= 1)
{
  if(!all(sapply(counts.list, is.data.table)))
    counts.list <- lapply(counts.list, as.data.table)
  if(is.null(names(counts.list)))
    names(counts.list) <- seq(counts.list)
  if(!is.null(control.cl) && any(!control.cl %in% names(counts.list)))
    stop("control.cl should match names(counts.list)")

  # Compute enrichment in non-control clusters
  cmb <- data.table(cl= names(counts.list))
  if(is.null(control.cl))
    cmb <- cmb[, .(ctl= cmb$cl[cmb$cl!=cl]), cl] else
      cmb <- cmb[, .(ctl= control.cl), cl]
  # Remove self-comparisons ----
  cmb[, self:= identical(cl, ctl), cl]
  cmb <- cmb[!(self), !"self"]
  # Motif enrichment ----
  res <- cmb[, {
    vl_motifEnrich(counts = counts.list[[cl]],
                   control.counts = rbindlist(counts.list[ctl]),
                   names= names,
                   plot= F)
  }, cl]
  res[, cl:= factor(cl, unique(cl))]
  setattr(res, "class", c("vl_enr_cl", "data.table", "data.frame"))

  # plot
  if(plot)
  {
    DT <- plot.vl_enr_cl(obj = res,
                         padj.cutoff= padj.cutoff,
                         top.enrich= top.enrich,
                         min.counts= min.counts,
                         order= order,
                         x.breaks= x.breaks,
                         color.breaks= color.breaks,
                         cex.balloons= cex.balloons,
                         col= col,
                         main= main,
                         plot.empty.clusters = plot.empty.clusters)
    if(add.motifs)
      vl_add_motifs(DT,
                    cex.width= cex.width,
                    cex.height= cex.height)
  }

  return(res)
}
