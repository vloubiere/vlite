heatmap.get.clusters <- function(name,
                                 idx,
                                 annot,
                                 x,
                                 clusters,
                                 distance,
                                 method,
                                 cutree,
                                 kmeans.k,
                                 cluster.seed,
                                 annot.col,
                                 cluster.col,
                                 gap.width)
{
  # Initiate object ----
  obj <- data.table(name= rownames(x),
                    idx= seq(nrow(x)),
                    annot = annot)

  # Define clustere ----
  if(length(clusters) == nrow(x)) {
    # User defined
    obj[, cluster:= clusters]
  } else if(isTRUE(clusters)) {
    # Clustering
    set.seed(cluster.seed)
    # Using kmeans
    if(!is.na(kmeans.k)) {
      obj[, cluster:= kmeans(x, centers = kmeans.k)$cluster]
    } else {
      # Compute distances for hclust
      .d <- if(distance %in% c("pearson", "spearman")) {
        as.dist(1 - cor(t(x),
                        use= "pairwise.complete.obs",
                        method= distance))
      } else {
        dist(x, method = distance)
      }
      # Hierachical clustering
      hcl <- hclust(.d, method = method)
      obj[, order:= hcl$order]
      # Cutree
      if(!missing(cutree)) {
        obj[, cluster:= cutree(hcl, cutree)]
      }
      # Extract dendrograms
      dend <- ggdendro::dendro_data(hcl,
                                    type = "rectangle",
                                    rotate= T)
      dend <- data.table::as.data.table(dend$segments)
    }
  } else if(!isFALSE(clusters))
    stop("Row and col clusters should match the dimensions of x or be logical vectors of length 1.")

  # Default values
  if(!"cluster" %in% names(obj))
    obj[, cluster:= as.character(NA)]
  obj[, cluster:= factor(cluster, sort(unique(cluster)))]
  if(!"order" %in% names(obj))
    obj[, order:= order(cluster)]

  # Add colors
  obj[, cluster.col:= colorRampPalette(cluster.col)(.NGRP)[as.numeric(cluster)], cluster]
  if(!is.null(annot))
    obj[, annot.col:= colorRampPalette(annot.col)(.NGRP)[as.numeric(annot)], annot]

  # Compute image coord and plotting positions
  obj <- obj[(order)]
  obj[, im:= .I]
  obj[, pos:= .I+(.GRP-1)*gap.width, cluster]

  # Return dendrogram
  return(list(obj= obj[order(idx)],
              dend= if(exists("dend")) dend else NULL))
}
