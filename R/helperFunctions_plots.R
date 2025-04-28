heatmap.get.clusters <- function(obj,
                                 x,
                                 clusters,
                                 distance,
                                 method,
                                 cutree,
                                 kmeans.k,
                                 cluster.seed,
                                 gap.width)
{
  # User defined clusters
  if(length(clusters) == nrow(x)) {
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
    stop("Row and col clusters should either match the dimensions of x or be logical vectors of length 1.")

  # Default values
  if(!"cluster" %in% names(obj))
    obj[, cluster:= as.character(NA)]
  if(!"order" %in% names(obj))
    obj[, order:= order(cluster)]

  # Compute plotting parameters
  obj <- obj[(order)]
  obj[, cluster:= factor(cluster, sort(unique(cluster)))]
  obj[, im:= .I]
  obj[, pos:= .I+(.GRP-1)*gap.width, cluster]

  # Return dendrogram
  return(list(obj= obj,
              dend= if(exists("dend")) dend else NULL))
}
