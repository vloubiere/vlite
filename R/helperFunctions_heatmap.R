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
  list(x= x, checkClass= checkClass, allLvls= if(checkClass=="factor") allLvls else NULL)
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
                                 som.grid,
                                 order.cl,
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
    
    if(!is.na(kmeans.k)) {
      
      # K-means clustering ----
      km <- kmeans(x, centers = kmeans.k)
      obj[, cluster:= factor(km$cluster)]
      obj[, order:= order(cluster)]
      
    } else if(!is.null(som.grid)) {
      
      # Self-organizing maps ----
      som.model <- vlite::somClustering(layers = x, grid = som.grid, init.seed = cluster.seed)
      obj[, cluster:= factor(som.model$unit.classif)]
      obj[, order:= order(cluster)]
      
    } else {
      
      # Hierachical clustering ----
      
      # Compute distances
      .d <- if(distance %in% c("pearson", "spearman")) {
        as.dist(1 - cor(t(x),
                        use= "pairwise.complete.obs",
                        method= distance))
      } else {
        dist(x, method = distance)
      }
      
      # Cluster
      hcl <- hclust(.d, method = method)
      obj[, order:= hcl$order]
      
      # Cutree
      if(!is.null(cutree)) {
        obj[, cluster:= cutree(hcl, cutree)]
        obj[, cluster:= factor(cluster)]
      }
      
      # Extract dendrograms ----
      dend <- ggdendro::dendro_data(hcl,
                                    type = "rectangle",
                                    rotate= T)
      dend <- data.table::as.data.table(dend$segments)
    }
  }
  
  # Order row clusters ----
  if(!is.null(order.cl)) {
    
    # Checks
    stopifnot(dim=="row")
    # Extract centers
    centers <- rbindlist(split(as.data.table(order.cl), obj$cluster), idcol = "cluster")
    centers <- centers[, lapply(.SD, mean), cluster]
    centers <- as.matrix(centers, 1)
    
    # zscore
    centers <- t(scale(t(centers), center= F))
    max.col <- apply(centers, 1, which.max)
    max.var <- apply(centers, 1, max)
    ordered <- levels(obj$cluster)[order(max.col, -max.var)]
    # Save new order
    obj[, cluster:= factor(cluster, ordered)]
    obj[, order:= order(cluster)]
  }
  
  # Add clusters colors ----
  if(length(cluster.col) < nlevels(obj$cluster))
    cluster.col <- colorRampPalette(cluster.col)(nlevels(obj$cluster))
  obj[, cluster.col:= cluster.col[cluster]]
  
  # Add annotations colors ----
  if(!is.null(annot)) {
    if(length(annot.col) < nlevels(obj$annot))
      annot.col <- colorRampPalette(annot.col)(nlevels(obj$annot))
    obj[, annot.col:= annot.col[annot]]
  }
  
  # Compute plotting coordinates ----
  obj <- obj[(order)]
  obj[, pos:= .I+(.GRP-1)*gap.width, cluster]
  
  # Compute dendrograms' plotting positions ----
  if(exists("dend") && is.null(order.cl)) {
      
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
