#' @describeIn vl_STRING_interaction Method to plot STRING interaction igraphs
#' @export
plot.vl_STRING <- function(obj, ...)
{
  .i <- stringToIgraph(obj)
  plot(.i, ...)
}
