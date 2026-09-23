#' Title
#'
#' @param x 
#' @param y 
#' @param ... 
#'
#' @returns
#' @export
#'
#' @examples
vl_fisher <- function(x, y= NULL, ...) {
  # Checks
  if(is.matrix(x)) {
    stopifnot(identical(dim(x), c(2,2)))
  } else {
    x <- table(x, y)
  }
  # Compute test
  .f <- fisher.test(x, y= NULL, ...)

  # Add corrected log2OR
  .f$log2OR.corr <- if(any(x==0)) {
    x <- x+.5
    log2((x[1,1] * x[2,2]) / (x[2,1] * x[1,2]))
  } else {
    log2(.f$estimate)
  }
  
  # Return
  return(.f)
}