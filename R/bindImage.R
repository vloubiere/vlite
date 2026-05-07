#' Title
#'
#' @param list 
#' @param dim 
#' @param frame 
#'
#' @returns
#' @export
#'
#' @examples
bindImage <- function(list, dim = 1) {
  stopifnot(dim %in% c(1, 2))
  
  # Check the dimension that must match
  check <- if(dim == 1)
    sapply(list, function(x) dim(x)[2]) else
      sapply(list, function(x) dim(x)[1])
  
  # Add lines/columns when missing to form a rectangle
  missing <- which(check < max(check))
  for(i in missing) {
    add <- max(check) - check[i]
    d <- dim(list[[i]])
    
    if(length(d) == 2) {
      list[[i]] <- if(dim == 1)
        cbind(list[[i]], matrix(0, nrow = d[1], ncol = add)) else
          rbind(list[[i]], matrix(0, nrow = add, ncol = d[2]))
    } else if(length(d) == 3) {
      pad <- if(dim == 1)
        array(0, dim = c(d[1], add, d[3])) else
          array(0, dim = c(add, d[2], d[3]))
      list[[i]] <- if(dim == 1)
        abind::abind(list[[i]], pad, along = 2) else
          abind::abind(list[[i]], pad, along = 1)
    }
  }
  
  # Do the binding
  tile_arr <- abind::abind(list, along = dim)
  tile <- EBImage::Image(tile_arr, colormode = "Color")
  
  # Return 
  return(tile)
}