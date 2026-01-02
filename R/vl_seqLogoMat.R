#' Plot logo from a numeric matrix
#'
#' Generic function that plots a logo out of an arbitrary numeric matrix. 
#'
#' @param mat A matrix (or a TFBSTools PFMatrix, PWMatrix or ICMatrix) with rows 'A', 'C', 'G', 'T'.
#' @param name The name of the motif which will be added as a title.
#' @param xlab xlab. Default= "Index".
#' @param ylab ylab. Default= "Frequency".
#' @param xlim xlim for the plot. Default= NULL.
#' @param ylim ylim for the plot. Default= NULL.
#' @param ... Extra parameters for the plot function.
#'
#' @examples
#' # Select br PPM
#' pfm.file <- system.file("extdata/hand_curated_Dmel_motifs_SCENIC_lite_Dec_2025.pfm", package = "vlite")
#' PPM <- importJASPAR(pfm.file, pseudocount = 0)$PPM[["br"]]
#'
#' # Plot
#' vl_par()
#' vl_seqLogoMat(PPM, name= "br")
#' 
#' @export
vl_seqLogoMat <- function(mat,
                          name= NULL,
                          xlab= "Index",
                          ylab= "Frequency",
                          xlim= NULL,
                          ylim= NULL,
                          ...)
{
  # Checks ----
  if(inherits(mat, c("PFMatrix", "PWMatrix", "ICMatrix")))
    mat <- mat@profileMatrix
  stopifnot(is.matrix(mat))
  stopifnot(identical(rownames(mat), c('A', 'C', 'G', 'T')))
  if(any(abs(mat)==Inf))
    stop("Infinite values not allowed. Use a pseudocount.")
  
  # Compute plotting position ----
  dat <- melt(as.data.table(mat, keep.rownames = T), id.vars = 'rn')
  setorderv(dat, c('variable', 'value'))
  dat[, left:= .GRP-1+.5, variable]
  dat[, ytop:= cumsum(value), variable]
  
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
      width = 1
    )
  }, .(rn, left, ytop, value)]
  
  # Add title ----
  if(!is.null(name))
    title(main= name)
}
