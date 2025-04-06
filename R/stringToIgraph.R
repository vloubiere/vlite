#' @describeIn vl_STRING_interaction Transforms a vl_STRING object into regular igraph
#' @export
stringToIgraph <- function(obj, directed= F)
{
  .i <- igraph::graph_from_data_frame(d = obj$E,
                                      vertices = obj$V,
                                      directed = directed)
  return(.i)
}
