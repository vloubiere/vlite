# Helper function to coerce each element of a list to numeric matrices
coerceLayerToNumeric <- function(x) {
  # Convert
  x <- as.data.table(x)
  
  # Select columns that can be contrasted
  sel.contrast <- which(sapply(x, function(x) length(unique(x))>1))
  x <- x[, sel.contrast, with= FALSE]
  
  # Factors/characters to one-hot-encoded
  res <- lapply(seq_len(ncol(x)), function(i) {
    y <- x[[i]]
    y <- factor(y)
    if (anyNA(y))
      y <- addNA(y)
    mm <- model.matrix(~ y - 1)
    colnames(mm) <- paste0(names(x)[i], "__", levels(y))
    mm
  })
  
  # Cbind
  final <- do.call(cbind, res)
  storage.mode(final) <- "numeric"
  
  # Return
  return(final)
}
