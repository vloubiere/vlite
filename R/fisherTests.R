#' fisherTests
#' 
#' Perform pairwise fisher tests between the elements of two lists.
#'
#' @param x.list A list of elements to be compared to y.list in the search of overlaps.
#' @param y.list A list of elements to be compared to x.list in the search of overlaps. 
#' @param universe If provided, will define the universe (background) of the test. Default
#' is to use the intersection between the (unique) values of x.list and y.list.
#' @param alternative The alternative to use for the fisher test. Default= "two.sided".
#'
#' @returns
#' @export
#'
#' @examples
fisherTests <- function(
    x.list,
    y.list,
    universe= NULL,
    alternative= "two.sided"
) {
  # Checks ----
  if(!is.list(x.list) && is.vector(x.list))
    x.list <- list(x= x.list)
  if(!is.list(y.list) && is.vector(y.list))
    y.list <- list(y= y.list)
  stopifnot(is.list(x.list))
  stopifnot(is.list(y.list))
  
  # Check names ----
  if(is.null(names(x.list)))
    names(x.list) <- paste0(seq_along(x.list), ".x")
  if(is.null(names(y.list)))
    names(y.list) <- paste0(seq_along(y.list), ".y")
  names(x.list) <- make.unique(names(x.list))
  names(y.list) <- make.unique(names(y.list))
  
  # Simplify ----
  x.list <- lapply(x.list, function(x) unique(unlist(x)))
  y.list <- lapply(y.list, function(x) unique(unlist(x)))
  
  # Retrieve universe ----
  if(is.null(universe)) 
    universe <- unique(intersect(unlist(x.list), unlist(y.list)))
  universe <- unique(universe)
  
  # Overlap lists with universe ----
  before <- lengths(c(x.list, y.list))
  x.list <- lapply(x.list, intersect, universe)
  y.list <- lapply(y.list, intersect, universe)
  after <- lengths(c(x.list, y.list))
  rm <- before-after
  if(any(rm>0))
    message(
      c(
        "Values missing from universe were removed:\n",
        paste0(c(names(x.list), names(y.list)), " -> ", rm, "/", before, " (", round((rm)/before*100), "%)", "\n")
      )
    )
    
  # Make comp table ----
  comp <- CJ(x= names(x.list), y= names(y.list))
  comp[, x:= factor(x, names(x.list))]
  comp[, y:= factor(y, names(y.list))]
  comp[, var.x:= .(list(x.list[[as.character(x)]])), x]
  comp[, var.y:= .(list(y.list[[as.character(y)]])), y]
  comp[, n.x:= lengths(var.x)]
  comp[, n.y:= lengths(var.y)]
  comp[, n.universe:= length(universe)]
  if(any(comp$n.x==comp$n.universe | comp$n.y==comp$n.universe))
    stop(
      "Invalid Fisher setup: at least one tested set equals the universe after filtering. ",
      "This often means `universe` was inferred from a selected gene set. ",
      "Provide the true background universe."
    )
  
  # Compute names with numbers ----
  comp[, name.x:= paste0(x, " (n= ", length(var.x[[1]]), ")"), x]
  setorderv(comp, "x")
  comp[, name.x:= factor(name.x, unique(name.x))]
  comp[, name.y:= paste0(y, " (n= ", length(var.y[[1]]), ")"), y]
  setorderv(comp, "y")
  comp[, name.y:= factor(name.y, unique(name.y))]
  
  # Compute enrichment ----
  comp[, c("log2OR.corr", "estimate", "p.value", "n.intersect", "intersect"):= {
    intersect <- intersect(var.x[[1]], var.y[[1]])
    c(
      vl_fisher(
        factor(universe %in% var.x[[1]], c(FALSE, TRUE)),
        factor(universe %in% var.y[[1]], c(FALSE, TRUE)),
        alternative= alternative
      )[c("log2OR.corr", "estimate", "p.value")],
      list(
        n.intersect= length(intersect),
        intersect= list(intersect)
      )
    )
  }, .(x, y)]
  # Compute log2OR
  comp[, log2OR:= log2(estimate)]
  # Adjust p-value
  comp[, padj:= p.adjust(p.value, "fdr")]
  
  # Return
  res <- comp[, .(name.x, name.y, n.intersect, n.x, n.y, n.universe, log2OR.corr, log2OR, padj, intersect)]
  setorderv(res, "name.x")
  return(res)
}