#' somClustering
#'
#' Wrapper around kohonen::supersom for quick SOM clustering of matrices..
#'
#' @param layers Either a matrix or a list of matrices to cluster.
#' @param grid The size of the grid. Default= c(2, 3).
#' @param reorder.cl Vector of cluster names used for reordering. Default= NULL.
#' @param clip.perc An integer vector of length 2 (or a list of length(layers)) specifying the upper and lower 
#' percentile (between [0-1]) used for the clipping of numeric layers prior to scaling and clustering. Ignored for
#' character/factor categorical layers. Default= NULL, meaning no clipping is performed.
#' @param clip.per.column If set to FALSE, clipping will be performed matrix-wide. Default= TRUE.
#' @param scale.layers A boolean value (or a vector of length(layers)) specifying whether numeric layers should be
#' scaled (z-scored) per column after clipping and prior to clustering. Setting it to TRUE will emphasizes patterns
#' between conditions/columns, while keeping it to FALSE (default) emphasizes differences in magnitude. Ignored for
#' character/factor categorical layers. Default= FALSE.
#' @param user.weights User-defined weights for the different layers. Length should be match the length of layers' list.
#' Default= 1 for all layers.
#' @param normalizeDataLayers Boolean indicating whether distance.weights should be calculated, in order to avoid some
#' layers to overwhelm others, simply because of the scale of the data points. If set to TRUE, these distance.weights
#' are applied first, and user weights are applied on top. If set to FALSE, user weights are applied to the data immediately.
#' Default= TRUE.
#' @param init.seed Seed used for random initialization.
#' @param maxNA.fraction The maxNA fraction tolerated for a given line. Default= 0L.
#'
#' @return
#' @export
#'
#' @examples
somClustering <- function(
    layers,
    grid= c(2, 3),
    reorder.cl= NULL,
    clip.perc= NULL,
    clip.per.column= TRUE,
    scale.layers= FALSE,
    user.weights= NULL,
    normalizeDataLayers= TRUE,
    init.seed= 1,
    maxNA.fraction= 0L
) {
  # Checks ----
  if(is.vector(layers) && is.numeric(layers))
    layers <- matrix(layers)
  if(is.data.table(layers) || is.matrix(layers) || is.data.frame(layers))
    layers <- list(layers)
  stopifnot(length(unique(sapply(layers, nrow)))==1)
  if(is.null(user.weights))
    user.weights <- rep(1, length(layers))
  stopifnot(length(user.weights)==length(layers))
  if(!is.list(clip.perc))
    clip.perc <- lapply(seq_along(layers), function(x) clip.perc)
  if(length(clip.perc)==1 & length(layers)>1)
    clip.perc <- lapply(seq_along(layers), function(x) clip.perc[[1]])
  stopifnot(length(clip.perc)==length(layers))
  if(length(scale.layers)==1 & length(layers)>1)
    scale.layers <- rep(scale.layers, length(layers))
  stopifnot(length(scale.layers)==length(layers))
  
  # Layers pre-processing and formatting ----
  layers <- lapply(seq_along(layers), function(i) {
    
    # Current layer
    var <- layers[[i]]
    
    # If current layer is numeric
    if(all(apply(var, 2, is.numeric))) {
      
      # Coerce to matrix
      var <- as.matrix(var)
      
      # Clip extreme values
      if(!is.null(clip.perc[[i]])) {
        if(clip.per.column) {
          
          # Per-column clipping
          var <- apply(var, 2, function(x) {
            lim <- quantile(x, clip.perc[[i]], na.rm= T)
            x[x<lim[1]] <- lim[1]
            x[x>lim[2]] <- lim[2]
            return(x)
          }
          )
        } else {
          
          # Matrix-wide clipping
          clip <- quantile(var, clip.perc[[i]], na.rm= T)
          var[var<clip[1]] <- clip[1]
          var[var>clip[2]] <- clip[2]
        }
      }
      
      # Scale
      if(scale.layers[i])
        var <- scale(var)
      
    } else {
      
      # If current layer contains characters/factors
      var <- coerceLayerToNumeric(var)
      # No clipping
      # No scaling
    }
    return(var)
  }
  )
  
  # Check layers formatting ----
  stopifnot(all(sapply(layers, is.numeric)))
  
  # Kohonen grid ----
  grid <- kohonen::somgrid(
    grid[1],
    grid[2],
    "hexagonal",
    toroidal= T
  )
  
  # Initialize grid ----
  uninf <- "first.attempt"
  while(uninf=="first.attempt" || isTRUE(uninf)) {
    # If uninformative layers were samples
    if(uninf != "first.attempt") {
      init.seed <- sample(1e3, size = 1)
      message(paste0("Trying to avoid uninformative layers with new seed ", init.seed))
    }
    # Layers initialization
    init <- lapply(layers, function(x)
    {
      set.seed(init.seed)
      x <- x[sample(nrow(x), grid$xdim*grid$ydim), , drop= F]
      return(x)
    })
    # Update number of unique combinations
    uninf <- any(sapply(init, function(x) nrow(unique(as.data.table(x))))==1)
  }
  
  # Clustering ----
  som <- kohonen::supersom(
    data = layers,
    grid= grid,
    init = init,
    user.weights= user.weights,
    maxNA.fraction = maxNA.fraction,
    normalizeDataLayers = normalizeDataLayers
  )
  
  # Reorder cl ----
  if(!is.null(reorder.cl)) {
    som$unit.classif <- factor(
      som$unit.classif,
      levels = c(reorder.cl, setdiff(unique(som$unit.classif), reorder.cl))
    )
  }
  
  # Return som object ----
  som$clip.perc <- clip.perc
  return(som)
}
