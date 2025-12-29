#' Plot PWM seqlogo
#'
#' Plots an ICM after removing low information bases (see ICM schneider correction).
#'
#' @param ICM A TFBSTools ICMatrix or ICMatrixList.
#' @param sel If ICM is a list, a vector of indices or of motif IDs specifying which motifs should be plotted. Default= 1.
#' @param name Name of the motif to plot as a title. By default, the name stored in the ICM will be used.
#'
#' @examples
#' # Retrieve br motif
#' pfm.file <- system.file("extdata/hand_curated_Dmel_motifs_SCENIC_lite_Dec_2025.pfm", package = "vlite")
#' mot <- importJASPAR(pfm.file)
#' idx <- which(mot$metadata$name=="br")
#' ICM <- mot$ICM[[idx]]
#'
#' # Plot
#' vl_seqLogo(ICM)
#'
#' @export
vl_seqLogo <- function(ICM,
                       sel= 1,
                       name= NULL)
{
  # Checks ----
  if(class(ICM)=="ICMatrixList") {
    if(is.integer())
  }

  if(is.null(name))
    name <- TFBSTools::name(ICM)

  # ICM to mat ----
  mat <- ICM@profileMatrix

  # Remove low information bases ----
  mat[, colSums(mat) < 0] <- 0
  mat[mat < 0] <- 0
  sel <- range(which(colSums(mat)>0))
  mat <- mat[, sel[1]:sel[2]]

  # Compute plotting positions ----
  dat <- as.data.table(mat, keep.rownames = "base")
  .c <- melt(dat, "base", variable.name = "xleft", value.name = "height")
  .c[, xleft:= as.numeric(xleft)]
  setorderv(.c, "height", 1)
  .c[, ytop:= cumsum(height), xleft]

  # Plot ----
  plot(NA,
       xlim= c(0.5, ncol(dat)+.5),
       ylim= c(0, max(.c$ytop)),
       xlab= "Position",
       ylab= "Bits",
       main= name)
  .c[, {
    plotDNAletter(base[1],
                  xleft[1],
                  ytop[1],
                  1,
                  height[1])
  }, .(base, xleft, height, ytop)]
  title(main= name)
}
