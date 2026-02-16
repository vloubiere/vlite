#' Clean pie chart
#'
#' @param x Variable to be summarized
#' @param col colors to be used
#' @param labels Type of labels to be added. "r"= number, "p"= percentage
#' @param ... extra arguments passed to pie(x, ...)
#'
#' @return Plots a pie chart
#' @export
#'
#' @examples
#' vl_pie(rep(1:3, 1:3))
vl_pie <- function(x,
                   col= NULL,
                   labels= "r",
                   remove.zeros= TRUE,
                   ...)
{
  .c <- table(x)
  if(is.null(col))
    col <- adjustcolor(rainbow(length(.c)), .3)
  if(labels=="r")
    names(.c) <- paste0(names(.c), " (", formatC(.c, big.mark = ","), ")")
  if(labels=="p")
    names(.c) <- paste0(names(.c), " (", round(.c/sum(.c)*100, 1), "%)")
  if(remove.zeros) {
    col <- col[.c>0]
    .c <- .c[.c>0]
  }
  pie(
    .c,
    col= col,
    ...
  )
}
