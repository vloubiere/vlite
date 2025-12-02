#' somClustering
#'
#' A wrapper around the som kohonen function to allow quick tests.
#'
#' @param layers Either a matrix or a list of matrices to cluster.
#' @param grid THe size of the grid. Default= c(2,3).
#' @param clip.perc Percentile used to clip numeric matrices. Set to c(0,1) to skip clipping.
#' @param scale Should the data be scaled using z-score? Default= FALSE
#' @param user.weights User-defined weights for the different layers. Length should be match the length of layers' list. Default= 1 for all layers.
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
    clip.perc= c(0.05, 0.95),
    scale= FALSE,
    user.weights= NULL,
    init.seed= 1,
    maxNA.fraction= 0L
) {
  # Checks
  if(is.numeric(layers))
    layers <- matrix(layers)
  if(is.data.table(layers) || is.matrix(layers) || is.data.frame(layers))
    layers <- list(layers)
  numeric.var <- sapply(layers, function(x) all(apply(x, 2, is.numeric)))
  layers <- coerce_layers_to_numeric(layers)
  stopifnot(length(unique(sapply(layers, nrow)))==1)
  stopifnot(all(sapply(layers, is.numeric)))
  if(is.null(user.weights))
    user.weights <- rep(1, length(layers))
  stopifnot(length(user.weights)==length(layers))

  # Scale ----
  if(scale) {
    layers <- layers[numeric.var] <- lapply(
      layers[numeric.var],
      function(x) {
        apply(
          x,
          2,
          scale
        )
      }
    )
  }

  # Clip outliers ----
  layers[numeric.var] <- lapply(
    layers[numeric.var],
    function(x) {
      apply(
        x,
        2,
        function(y) {
          lim <- quantile(y, clip.perc, na.rm= T)
          y[y<lim[1]] <- lim[1]
          y[y>lim[2]] <- lim[2]
          return(y)
        }
      )
    }
  )

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
    maxNA.fraction = maxNA.fraction
  )

  # Return som object ----
  return(som)
}
