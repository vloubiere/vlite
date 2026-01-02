#' Add ICM logos to an existing plot
#'
#' plot seqlogo from pwm matrix
#'
#' @param ICM A list of ICMs to be plotted.
#' @param sel Indices or names of the ICM that should be plotted. Should match the length of x and y.
#' Default= NULL.
#' @param x x positions (see pos argument).
#' @param y y positions (centered).
#' @param pos Either 2 (left) of 4 (right)
#' @param cex.width width expansion factor.
#' @param cex.height height expansion factor.
#' @param min.content Flanks with a smaller summed content are not plotted.
#'
#' @examples
#' # Retrieve br motif
#' pfm.file <- system.file("extdata/hand_curated_Dmel_motifs_SCENIC_lite_Dec_2025.pfm", package = "vlite")
#' ICM <- importJASPAR(pfm.file)$ICM
#'
#' # Plot
#' plot(0, 1, type= "n")
#' addSeqLogo(
#' ICM, 
#' x = -1,
#' y= seq(0.8, 1.2, .2),
#' sel= c("br", "dl", "Eip74EF"),
#' cex.height = 1.5,
#' cex.width = 1, 
#' pos= 2
#' )
#'
#' @export
addSeqLogo <- function(ICM,
                       sel= NULL,
                       x= NULL,
                       y= NULL,
                       pos= 2,
                       cex.width= 1,
                       cex.height= 1)
{
  # Checks ----
  if(class(ICM)=="ICMatrix")
    ICM <- TFBSTools::ICMatrixList(ICM)
  stopifnot(inherits(ICM, "ICMatrixList"))
  if(!is.null(sel))
    ICM <- ICM[sel]
  if(length(x)==1)
    x <- rep(x, length(ICM))
  if(length(y)==1)
    y <- rep(y, length(ICM))
  stopifnot(length(x)==length(ICM))
  stopifnot(length(y)==length(ICM))
  
  # Compute width and height ----
  .w <- strwidth("M", cex= cex.width)
  .h <- strheight("M", cex= cex.height)
  
  # Loop ----
  lapply(seq_along(ICM), function(i) {
    # ICM to mat ----
    mat <- ICM[[i]]@profileMatrix
    
    # Remove low information bases ----
    mat[, colSums(mat) < 0] <- 0
    mat[mat < 0] <- 0
    sel <- range(which(colSums(mat)>0))
    mat <- mat[, sel[1]:sel[2]]
    
    # Compute plotting positions ----
    dat <- as.data.table(mat, keep.rownames = "base")
    .c <- melt(dat, "base", variable.name = "xleft", value.name = "height")
    .c[, xleft:= x[i]+as.numeric(xleft)*.w]
    setorderv(.c, "height", 1)
    .c[, height:= height/max(height)*.h]
    .c[, ytop:= y[i]+cumsum(height)-.h/2, xleft]
    
    # Revert ----
    if(pos==4)
      .c[, xleft:= xleft-(ncol(mat)+2)*.w]
    
    # Plot ----
    .c[, {
      plotDNAletter(base[1],
                    xleft[1],
                    ytop[1],
                    .w,
                    height[1])
    }, .(base, xleft, height, ytop)]
  }
  )
}
