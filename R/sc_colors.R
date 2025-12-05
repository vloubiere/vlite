#' Plot and cluster top marker genes from single-cell transcriptome data
#'
#' Default colors used for cell cluster in sc UMAP
#'
#' @param clusters A factor vector containin clusters to color.
#' @param unique.lvl If set to TRUE, one color is returned per unique level inside clusters.
#'
#' @return The default colors of ggplot
#'
#' @export
sc_colors <- function(clusters, unique.lvl= FALSE)
{
  # Checks ----
  if(!is.factor(clusters))
    stop("clusters should be a factor vector.")
  if(unique.lvl)
    clusters <- factor(levels(clusters), levels(clusters))

  # Return ----
  return(scales::hue_pal()(length(levels(clusters)))[clusters])
}
