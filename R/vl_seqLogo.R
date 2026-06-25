#' Plot motif logo from a PFM, PPM or ICM matrix
#'
#' Generic function that plots a motif logo out of a matrix. 
#'
#' @param mat A matrix (or a TFBSTools PFMatrix, PWMatrix or ICMatrix) with rows 'A', 'C', 'G', 'T'.
#' @param name The name of the motif which will be added as a title.
#' @param xlab xlab. Default= "Position".
#' @param ylab ylab. Default= "Frequency" for PFMatrix, "Probability" for PWMatrix and "Bits" for ICMatrix.
#' @param xlim xlim for the plot. Default= NULL.
#' @param ylim ylim for the plot. Default= NULL.
#' @param ... Extra parameters for the plot function.
#'
#' @examples
#' # Select br PPM
#' pfm.file <- system.file("extdata/hand_curated_Dmel_motifs_SCENIC_lite_Dec_2025.pfm", package = "vlite")
#' mot <- vlite::importJASPAR(pfm.file)
#'
#' # Plot
#' vl_par()
#' vl_seqLogo(mot$PFM[[1]])
#' vl_seqLogo(mot$PPM[[1]])
#' vl_seqLogo(mot$ICM[[1]])
#' 
#' @export
vl_seqLogo <- function(
    mat,
    name= NULL,
    xlab= "Position",
    ylab= NULL,
    xlim= NULL,
    ylim= NULL,
    ...
)
{
  # Checks ----
  if(inherits(mat, c("PFMatrix", "PWMatrix", "ICMatrix"))) {
    if(is.null(name))
      name <- mat@name
    if(is.null(ylab))
      ylab <- switch(
        class(mat)[1],
        "PFMatrix"= "Frequency",
        "PWMatrix"= "Probability",
        "ICMatrix"= "Bits"
      )
    mat <- mat@profileMatrix
  }
  if(min(mat, na.rm = T) < 0)
    stop("This method should not be used to represent matrices with negative values (e.g. log2probratio PWM...)")
  stopifnot(is.matrix(mat))
  stopifnot(identical(rownames(mat), c('A', 'C', 'G', 'T')))
  if(any(abs(mat)==Inf))
    stop("Infinite values not allowed. Use a pseudocount.")
  
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
  dat[, left:= .GRP-1+.5, variable]
  dat[, ytop:= cumsum(value), variable]
  dat[, width:= 1, variable]
  
  # Compute xlim and ylim ----
  if(is.null(xlim))
    xlim <- c(0.5, ncol(mat)+.5)
  if(is.null(ylim))
    ylim <- c(0, max(dat$ytop, na.rm = T))
  
  # Initiate plot ----
  plot(
    NA,
    type= 'n',
    xlim= xlim,
    ylim= ylim,
    xlab= xlab,
    ylab= ylab,
    ...
  )
  
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
  
  # Add title ----
  if(!is.null(name))
    title(main= name)
}
