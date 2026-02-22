#' Title
#'
#' @param x 
#' @param col 
#' @param labels labels to be added to the lines (right part of the plot).
#' @param names The condition names to be plotted below the x axis. Default= colnames(x).
#' @param xlim 
#' @param ylim 
#' @param ... 
#'
#' @returns
#' @export
#'
#' @examples
profilePlot <- function(
    x,
    col,
    labels= NULL,
    names= colnames(x),
    ylab= "Variable",
    xlim= NULL,
    ylim= NULL,
    main= NA,
    ...
) {
  if(!is.matrix(x))
    x <- as.matrix(x)
  if(!is.numeric(x))
    stop("x should be a numeric matrix")
  if(is.function(names))
    names <- sapply(colnames(x), function(x) names(x))
  
  # Initiate plot
  plot(
    c(1, ncol(x)),
    range(c(x), na.rm= T),
    type= "n",
    xlim= xlim,
    ylim= ylim,
    xlab= NA,
    ylab= ylab,
    xaxt= "n",
    main= main
  )
  vlite::tiltAxis(seq(ncol(x)), labels = names)
  
  # Plot lines
  sapply(seq(nrow(x)), function(i) lines(x[i,], col= col[i], ...))
  
  # Add labels
  text(ncol(x), x[, ncol(x)], labels, col= col, pos= 4, xpd= NA, cex= .7)
}