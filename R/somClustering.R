#' somClustering
#'
#' A wrapper around the som kohonen function to allow quick tests.
#'
#' @param layers Either a matrix or a list of matrices to cluster.
#' @param grid The size of the grid. Default= c(2,3).
#' @param clip.perc An integer vector of length 2 (or a list of such vectors) specifying the upper and lower 
#' percentile used to clip numeric layers data before clustering (ignored for characters/factors categorical layers).
#' Default= NULL, meaning no clipping is performed.
#' @param scale Should numeric layers be scaled (z-scored) before clustering? Default= FALSE.
#' @param user.weights User-defined weights for the different layers. Length should be match the length of layers' list. Default= 1 for all layers.
#' @param normalizeDataLayers Boolean indicating whether distance.weights should be calculated (see details section).
#' If set to FALSE, user weights are applied to the data immediately. Default= TRUE.
#' @param init.seed Seed used for random initialization/
#' @param maxNA.fraction The maxNA fraction tolerated for a given line. Default= 0L.
#'
#' @return
#' @export
#'
#' @examples
somClustering <- function(
    layers,
    grid= c(2, 3),
    clip.perc= c(0, 1),
    scale= FALSE,
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
  
  # Layers pre-processing and formatting ----
  layers <- lapply(seq_along(layers), function(i) {
    
    # Current layer
    var <- layers[[i]]
    
    # If current layer is numeric
    if(all(apply(var, 2, is.numeric))) {
      
      # Coerce to matrix
      var <- as.matrix(var)
      
      # Scale
      if(scale)
        var <- apply(var, 2, scale)
      
      # Clip
      if(!is.null(clip.perc[[i]]))
        var <- apply(var, 2, function(x) {
          lim <- quantile(x, clip.perc[[i]], na.rm= T)
          x[x<lim[1]] <- lim[1]
          x[x>lim[2]] <- lim[2]
          return(x)
        }
        )
      
    } else {
      
      # If current layer contains characters/factors
      var <- coerceLayerToNumeric(var)
      # No scaling
      # No clipping
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
  init <- lapply(layers, function(x)
  {
    set.seed(init.seed)
    x <- x[sample(nrow(x), grid$xdim*grid$ydim), , drop= F]
    return(x)
  })
  
  # Clustering ----
  som <- kohonen::supersom(
    data = layers,
    grid= grid,
    init = init,
    user.weights= user.weights,
    maxNA.fraction = maxNA.fraction,
    normalizeDataLayers = normalizeDataLayers
  )
  
  # Return som object ----
  som$clip.perc <- clip.perc
  return(som)
}
