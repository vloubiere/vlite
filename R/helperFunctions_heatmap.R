# Take rows and columns names/idx/annot and makes an object containing
# clusters, colors, positions...
heatmap.get.clusters <- function(dim= "row",
                                 x = x,
                                 annot,
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
  # Checks ----
  if(!dim %in% c("row", "col"))
    stop("dim should be one of row or col")
  if(dim=="col" & !is.na(kmeans.k))
    stop("Kmeans is only implemented for rows, not columns")
  if(dim=="col")
    x <- t(x)
  if(!is.null(annot)) {
    if(is.numeric(annot) && any(annot %% 1) != 0)
      stop("Row and col annotations can only contain integer or non-numeric values that will be coerced to factors.")
    annot <- factor(annot, sort(unique(annot)))
  }
  if(length(clusters) != nrow(x) && !(isFALSE(clusters) | isTRUE(clusters)))
    stop("Row and col clusters should match the dimensions of x or be logical vectors of length 1.")

  # Initiate object ----
  obj <- data.table(name= rownames(x),
                    idx= seq(nrow(x)),
                    annot = annot,
                    cluster = if(length(clusters) == nrow(x)) clusters else NA_character_,
                    order = seq(nrow(x)))

  # Clustering ----
  if(isTRUE(clusters)) {
    set.seed(cluster.seed)

    # Kmeans clustering
    if(!is.na(kmeans.k)) {
      obj[, cluster:= kmeans(x, centers = kmeans.k)$cluster]

    # Hierachical clustering
    } else {
      # Compute distances
      .d <- if(distance %in% c("pearson", "spearman")) {
        as.dist(1 - cor(t(x),
                        use= "pairwise.complete.obs",
                        method= distance))
      } else {
        dist(x, method = distance)
      }

      # Actual clustering
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
  }

  # Coerce to factors and compute order ----
  obj[, cluster:= factor(cluster, sort(unique(cluster)))]
  obj[, order:= order(cluster)]

  # Add clusters and annotations colors ----
  obj[, cluster.col:= colorRampPalette(cluster.col)(nlevels(cluster))[cluster]]
  if(!is.null(annot))
    obj[, annot.col:= colorRampPalette(annot.col)(nlevels(annot))[annot]]

  # Compute plotting coordinates ----
  obj <- obj[(order)]
  obj[, pos:= .I+(.GRP-1)*gap.width, cluster]

  # For dendrogram, compute plotting positions ----
  if(exists("dend")) {
    # interpolate cluster gaps
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
    # Normalized positions
    if(dim=="row") {
      dend[, x0:= y/diff(range(c(y, yend)))]
      dend[, x1:= yend/diff(range(c(y, yend)))]
      setnames(dend,
               c("s.start", "s.end"),
               c("y0", "y1"))
      dend[, y0:= y0-0.5]
      dend[, y1:= y1-0.5]
    } else {
      setnames(dend,
               c("s.start", "s.end"),
               c("x0", "x1"))
      dend[, y0:= y/diff(range(c(y, yend)))]
      dend[, y1:= yend/diff(range(c(y, yend)))]
    }
  } else {
    dend <- NULL
  }

  # Back to original order (don't do it earlier!) ----
  obj <- obj[order(idx)]

  # Set rownames ----
  if(dim=="row") {
    setnames(obj,
             c("idx", "pos"),
             c("line.idx", "y.pos"))
  } else {
    setnames(obj,
             c("idx", "pos"),
             c("column.idx", "x.pos"))
  }

  # Return object and dendrogram ----
  return(list(obj= obj,
              dend= dend))
}
