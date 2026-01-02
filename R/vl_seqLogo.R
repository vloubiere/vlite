#' Plot PWM seqlogo
#'
#' Plots an ICM after removing low information bases (see ICM schneider correction).
#'
#' @param ICM A TFBSTools ICMatrix or ICMatrixList.
#' @param sel If ICM is an ICMatrixList, index or name of the ICM that should be plotted.
#' Default= 1.
#' @param name Name of the motif to plot as a title. By default, the name stored in the ICM will be used.
#' @param ... Extra plotting parameters for the plot funciton.
#'
#' @examples
#' # Import motifs ICMs
#' pfm.file <- system.file("extdata/hand_curated_Dmel_motifs_SCENIC_lite_Dec_2025.pfm", package = "vlite")
#' ICM <- importJASPAR(pfm.file, pseudocount = 0)$ICM
#'
#' # Plot
#' vl_seqLogo(ICM, sel= "br")
#'
#' @export
vl_seqLogo <- function(ICM,
                       sel= 1,
                       name= NULL,
                       ...)
{
  # Checks ----
  if(inherits(ICM, "ICMatrixList")) {
    stopifnot(length(sel)==1)
    ICM <- ICM[[sel]]
  }
  stopifnot(inherits(ICM, "ICMatrix"))
  if(is.null(name))
    name <- TFBSTools::name(ICM)

  # ICM to mat ----
  mat <- ICM@profileMatrix

  # Remove low information bases ----
  mat[, colSums(mat) < 0] <- 0
  mat[mat < 0] <- 0
  sel <- range(which(colSums(mat)>0))
  mat <- mat[, sel[1]:sel[2]]

  # Plot ----
  vl_seqLogoMat(
    mat,
    name = name,
    xlab= "Position",
    ylab= "Bits"
  )
}
