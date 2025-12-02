# Helper function to coerce each element of a list to numeric matrices
coerce_layers_to_numeric <- function(layers) {
  layers <- lapply(
    layers,
    function(x) {
      # Matrix
      if(is.matrix(x)) {
        return(x)
      } else if(all(apply(x, 2, is.numeric))) {
      # Numeric columns to matrix
        return(as.matrix(x))
      } else {
        # Factors/characters to one-hot-encoded
        .c <- lapply(
          seq_len(ncol(x)),
          function(i) {
            y <- x[[i]]
            if (!is.factor(y))
              y <- factor(y)
            if (anyNA(y))
              y <- addNA(y)
            mm <- model.matrix(~ y - 1)
            colnames(mm) <- paste(colnames(x)[i], levels(y), sep = "__")
            return(mm)
          }
        )
        # Cbind
        .c <- do.call(cbind, .c)
        storage.mode(.c) <- "numeric"
        return(.c)
      }
    }
  )
  return(layers)
}
