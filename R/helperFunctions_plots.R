
# Take rows and columns names/idx/annot and makes an object containing
# clusters, colors, positions...
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

  # Define clusters ----
  if(isFALSE(clusters)) {
    # No clustering
    obj[, cluster:= as.character(NA)]
  } else if(length(clusters) == nrow(x)) {
    # User defined
    obj[, cluster:= clusters]
  } else if(isTRUE(clusters)) {
    # Clustering
    set.seed(cluster.seed)
    # Using kmeans
    if(!is.na(kmeans.k)) {
      obj[, cluster:= kmeans(x, centers = kmeans.k)$cluster]
    } else {
      # Distances
      .d <- if(distance %in% c("pearson", "spearman")) {
        as.dist(1 - cor(t(x),
                        use= "pairwise.complete.obs",
                        method= distance))
      } else {
        dist(x, method = distance)
      }
      # Hierachical clustering
      hcl <- hclust(.d, method = method)
      obj[, cluster:= as.character(NA)]
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
  } else
    stop("Row and col clusters should match the dimensions of x or be logical vectors of length 1.")

  # Coerce to factors and compute order ----
  obj[, cluster:= factor(cluster, unique(sort(cluster)))]
  obj[, cluster:= factor(cluster, unique(sort(cluster)))]
  if(!"order" %in% names(obj))
    obj[, order:= order(cluster)]

  # Compute image coordinates and plotting positions ----
  obj <- obj[(order)]
  obj[, im:= .I]
  obj[, pos:= .I+(.GRP-1)*gap.width, cluster]

  # Add clusters and annotations colors ----
  obj[, cluster.col:= colorRampPalette(cluster.col)(length(levels(cluster)))[as.numeric(cluster)]]
  if(!is.null(annot))
    obj[, annot.col:= colorRampPalette(annot.col)(length(levels(annot)))[as.numeric(annot)]]

  # For dendrogram, compute plotting positions by interpolating gaps
  dend <- if(exists("dend")) {
    dend[, s.start := {
      start <- obj$pos[floor(x)]
      end   <- obj$pos[ceiling(x)]
      start + (end - start) * (x - floor(x))
    }]
    dend[, s.end := {
      start <- obj$pos[floor(xend)]
      end   <- obj$pos[ceiling(xend)]
      start + (end - start) * (xend - floor(xend))
    }]
  } else
    NULL

  # Return object and dendrogram ----
  return(list(obj= obj[order(idx)], # Original order
              dend= dend))
}
