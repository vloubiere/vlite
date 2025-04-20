#' motifPosToMatrix
#'
#' Transforms the output of vl_motifPos() to a list of scores' matrices.
#'
#' @param motifPos A ?vl_motifPos() output data.table.
#' @param seqWidth. The width of the sequences for which motif positions have been computed.
#' @param cleanup.cache If set to TRUE, overwrites cached intermediate results. Default= FALSE.
#'
#' @return Returns a list of maximum scores' matrices (one per motif).
#' @export
#'
#' @examples
motifPosToMatrix <- function(motifPos,
                             seqWidth,
                             cleanup.cache= FALSE)
{
  # Checks
  if(!is.data.table(motifPos))
    stop("motifPos should be a data.table similar to ?vl_motifPos() output.")
  if(!all(c("motif", "seqlvls", "mot.count", "ir") %in% names(motifPos)))
    stop("motifPos should contain 'motif', 'seqlvls', 'mot.count' and 'ir' columns (see ?vl_motifPos()).")
  if(max(rbindlist(motifPos$ir)$end)>seqWidth)
    stop("Some coordinates within the motifPos object are higher than seqWidth.")

  # Use a temp directory for caching signal files ----
  cache_dir <- tempdir()
  params <- list(motifPos,
                 seqWidth)
  key <- digest::digest(params)
  output.file <- file.path(cache_dir, paste0(key, ".rds"))

  # Main function
  if(!file.exists(output.file) | cleanup.cache) {
    # Copy for encapsulation
    dat <- data.table::copy(motifPos)

    # Sort based on score (nts with several scores will receive the highest one)
    dat[!is.na(mot.count), ir:= lapply(ir, setorderv, "score")]

    # Initiate an empty vector
    empty <- rep(NA, seqWidth)

    # For each motif
    res <- dat[, {
      # For each seqlvls
      .c <- .SD[, {
        # Initiate an empty vector
        mat <- empty
        if(!is.na(mot.count)) {
          # Replace NAs with score
          ir[[1]][, {
            mat[start:end] <<- score
          }, .(start, end, score)]
        }
        .(.(mat))
      }, seqlvls]
      # Combine vectors into a unique matrix per motif
      mat <- do.call(rbind, .c$V1)
      rownames(mat) <- seqlvls
      # Return as a nested matrix
      .(mat= .(mat))
    }, motif]

    # Make as a list
    final <- res$mat
    names(final) <- res$motif

    # Save in cache
    saveRDS(final, output.file)
  } else
    final <- readRDS(output.file)

  # Return
  return(final)
}
