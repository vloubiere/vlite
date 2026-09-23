# Compute empirical cumulative recovery and AUC ----
cumRecov <- function(labels, ranks.matrix, max.rank, return.recov= T, decreasing= decreasing) {
  # For each rank
  res <- ranks.matrix[, lapply(.SD, function(x) {
    # Reorder labels based on ranks
    lab <- labels[order(x, decreasing = decreasing)][1:max.rank]
    # Compute recovery
    recov <- cumsum(lab)/sum(labels)
    AUC <- sum(recov)/length(recov)
    # Return
    if(return.recov)
      list(AUC= AUC, recov= recov) else
        AUC
  })]
  
  # If recovery curves should be returned
  if(return.recov) {
    # Aggregate result
    AUC <- unlist(res[1,])
    recov <- data.table::transpose(res[2,])$V1
    recov <- do.call(cbind, recov)
    colnames(recov) <- names(ranks.matrix)
    # Return
    final <- list(AUC= AUC, recov= recov)
    return(final)
  } else {
    # Simplify
    AUC <- unlist(res[1,])
    # Return
    return(AUC)
  }
}

# Return cancidate regions ----
candidateRegions <- function(labels, ranks.matrix, leading.edges, max.rank, decreasing= decreasing) {
  # Retrieve motifs with a leading edge
  res <- ranks.matrix[, !is.na(leading.edges), with= F]
  leading.edges <- na.omit(leading.edges)
  names(leading.edges) <- colnames(res)
  # Melt
  res <- melt(res, measure.vars= names(res), variable.name = "motif", value.name = "rank")
  # Select positively labelled ranks
  idx <- which(labels)
  res <- res[, c(.(row.idx= idx), .SD[idx]), motif]
  res <- res[, .SD[rank<=leading.edges[motif]], motif]
  setorderv(res, c("motif", "rank"), c(1, ifelse(decreasing, -1, 1)))
  # Return
  return(res)
}