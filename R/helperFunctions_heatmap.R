# Coerce x to numeric matrix (useful for factors) ----
toNumMatrix <-  function(x)
{
  if(!is.matrix(x)) {
    checkClass <- unique(sapply(x, class))
    if(length(checkClass)>1)
      stop("x cannot contain mixed classes.")
    if(checkClass=="factor") {
      allLvls <- unique(c(sapply(x, levels)))
      # Factors to numeric
      x <- lapply(x, function(x) as.numeric(factor(x, allLvls)))
      x <- do.call(cbind, lapply(x, as.numeric))
    }
    # Numeric matrix
    x <- as.matrix(x)
  } else {
    checkClass <- "non-factor"
  }
  if(is.logical(x))
    x <- apply(x, 2, as.numeric)
  if(!is.numeric(x))
    stop("x should contain numeric values, factors or logical values.")

  # Return matrix and checkClass
  list(x= x, checkClass= checkClass)
}

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
  if(dim=="col")
    x <- t(x)

  # Initiate object ----
  obj <- data.table(name= rownames(x),
                    idx= seq(nrow(x)),
                    annot = annot,
                    cluster = if(is.factor(clusters)) clusters else factor(NA),
                    order = seq(nrow(x)))
  obj[, order:= order(cluster)]

  # Clustering ----
  if(isTRUE(clusters)) {
    set.seed(cluster.seed)

    # Kmeans clustering
    if(!is.na(kmeans.k)) {
      obj[, cluster:= kmeans(x, centers = kmeans.k)$cluster]
      obj[, cluster:= factor(cluster)]
      obj[, order:= order(cluster)]

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
        obj[, cluster:= factor(cluster)]
      }

      # Extract dendrograms
      dend <- ggdendro::dendro_data(hcl,
                                    type = "rectangle",
                                    rotate= T)
      dend <- data.table::as.data.table(dend$segments)
    }
  }

  # Add clusters and annotations colors ----
  obj[, cluster.col:= colorRampPalette(cluster.col)(nlevels(cluster))[cluster]]
  if(!is.null(annot))
    obj[, annot.col:= colorRampPalette(annot.col)(nlevels(annot))[annot]]

  # Compute plotting coordinates ----
  obj <- obj[(order)]
  obj[, pos:= .I+(.GRP-1)*gap.width, cluster]

  # Compute dendrograms' plotting positions ----
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

  # Order based on position ----
  obj <- obj[order(pos)] # Don't do it earlier!

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
