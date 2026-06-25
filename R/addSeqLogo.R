#' Plot logo from a numeric matrix
#'
#' Generic function that plots a logo out of an arbitrary numeric matrix. 
#'
#' @param mat A matrix (or a TFBSTools PFMatrix, PWMatrix or ICMatrix) with rows 'A', 'C', 'G', 'T'.
#' @param x Th x position where the logo should start or finish (depending on 'pos' argument).
#' @param y The y position of the bottom part of the logo. Default= 0.
#' @param cex.width A width expansion factor applied to letter widths.
#' @param cex.width A width expansion factor applied to letter heights.
#' @param pos Specifies on which side of the x position the logo should be plotted. Wither 2 (left) or 4 (right). Default= 4.
#'
#' @examples
#' # Select br PPM
#' pfm.file <- system.file("extdata/hand_curated_Dmel_motifs_SCENIC_lite_Dec_2025.pfm", package = "vlite")
#' mot <- vlite::importJASPAR(pfm.file)
#'
#' # Plot
#' vl_par()
#' plot(0, 0, type= "n")
#' addSeqLogo(mot$PFM[[1]], x= -1, y= .5)
#' addSeqLogo(mot$PPM[[1]], x= -0.5, y= 0)
#' addSeqLogo(mot$ICM[[1]], x= 0, y= -0.5)
#' 
#' @export
addSeqLogo <- function(
    mat,
    x,
    y,
    cex.width= 1,
    cex.height= 1,
    pos= 4
)
{
  # Checks ----
  if(inherits(mat, c("PFMatrix", "PWMatrix", "ICMatrix"))) {
    mat <- mat@profileMatrix
  }
  if(min(mat, na.rm = T) < 0)
    stop("This method should not be used to represent matrices with negative values (e.g. log2probratio PWM...)")
  stopifnot(is.matrix(mat))
  stopifnot(identical(rownames(mat), c('A', 'C', 'G', 'T')))
  if(any(abs(mat)==Inf))
    stop("Infinite values not allowed. Use a pseudocount.")
  stopifnot(pos %in% c(2, 4))
  
  # Remove low information bases (ICM only) ----
  if(class(mat)[1]=="ICMatrix") {
    mat[, colSums(mat) < 0] <- 0
    mat[mat < 0] <- 0
    sel <- range(which(colSums(mat)>0))
    mat <- mat[, sel[1]:sel[2]]
  }
  
  # Compute plotting position ----
  dat <- melt(as.data.table(mat, keep.rownames = T), id.vars = 'rn')
  setorderv(dat, c('variable', 'value'))
  dat[, width:= strwidth("M")*cex.width]
  dat[, left:= (.GRP-1)*width+x, variable]
  max.height <- max(dat[, sum(value), variable]$V1)
  dat[, value:= value/max.height*(strheight("M")*(cex.height*2))]
  dat[, ytop:= cumsum(value)+y, variable]
  if(pos==2)
    dat[, left:= left-diff(range(left))-width]
  
  # Plot ----
  dat[, {
    vlite::plotDNAletter(
      letter = rn,
      xleft = left,
      ytop = ytop,
      height = value,
      width = width
    )
  }, .(rn, left, ytop, value, width)]
}
