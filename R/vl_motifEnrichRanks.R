#' Title
#'
#' @param labels A logical vector of labels matching the number of rows in ranks.matrix.
#' @param ranks.matrix A numeric matrix of motif ranks (one motif per column).
#' @param max.rank The maximum rank up to which cumulative recovery should be computed. Default= 5000.
#' @param N.permutations The number of permutations used for empirical p-values. 
#' If set to NA, p-values are not computed. Default= 1,000.
#' @param decreasing Should the provided ranks be used in decreasing order? Default = FALSE.
#' @param cleanup.cache Should the cache be cleaned? Default= FALSE.
#'
#' @returns
#' @export
#'
#' @examples
vl_motifEnrichRanks <- function(
    labels,
    ranks.matrix,
    max.rank= 5000,
    N.permutations= 1000,
    decreasing= F, 
    seed= 1,
    cleanup.cache = FALSE,
    tmp.dir= "/zssd/scratch/vincent.loubiere/projects/iCis_lite/db/tmp/"
) {
  # Checks ----
  if(is.matrix(ranks.matrix))
    ranks.matrix <- as.data.table(ranks.matrix)
  stopifnot(is.data.table(ranks.matrix))
  stopifnot(!is.null(names(ranks.matrix)))
  stopifnot(!any(duplicated(names(ranks.matrix))))
  stopifnot(nrow(ranks.matrix)==length(labels))
  stopifnot(nrow(ranks.matrix)>=max.rank)
  stopifnot(is.logical(labels))
  stopifnot(any(labels))
  stopifnot(!anyNA(labels))
  
  # Check if computed ----
  tmp <- vlite::vl_cache_file(
    input.list = list(
      labels,
      ranks.matrix,
      max.rank,
      N.permutations,
      decreasing,
      seed
    ),
    tmp.dir = tmp.dir
  )
  if(!file.exists(tmp) | cleanup.cache) {
    # Compute recovery ----
    res <- cumRecov(
      labels = labels,
      ranks.matrix = ranks.matrix,
      max.rank = max.rank,
      decreasing = decreasing
    )
    
    # Compute NES ----
    enr <- data.table(motif= names(ranks.matrix), AUC= res$AUC)
    enr[, NES:= (AUC - mean(AUC)) / sd(AUC)]
    
    # Leading edge ----
    mean.recov <- apply(res$recov, 1, mean)
    sd.recov <- apply(res$recov, 1, sd)
    enr$leading.edge <- apply(res$recov, 2, function(x) {
      delta <- x-(mean.recov+2*sd.recov)
      if(max(delta)>0) which.max(delta) else NA_integer_
    }
    )
    candidate_regions <- candidateRegions(
      labels = labels,
      leading.edges= enr$leading.edge,
      ranks.matrix = ranks.matrix,
      max.rank = max.rank,
      decreasing = decreasing
    )
    
    # Empirical p-value ----
    if(!is.na(N.permutations)) {
      set.seed(seed)
      rdm <- lapply(1:N.permutations, function(i) {
        cumRecov(
          labels = sample(labels, length(labels)),
          ranks.matrix = ranks.matrix,
          max.rank = max.rank,
          decreasing = decreasing,
          return.recov= FALSE
        )
      }
      )
      rdm <- do.call(rbind, rdm)
      enr$pval <- sapply(seq_along(res$AUC), function(i) (1+sum(rdm[,i]>=res$AUC[i]))/(nrow(rdm)+1))
      enr$padj <- p.adjust(enr$pval, "fdr")
    }
    
    # Format output ----
    output <- list(
      enr= enr,
      recov= list(
        mean.recov= mean.recov,
        sd.recov= sd.recov,
        motif.recov= res$recov
      ),
      candidate_regions= candidate_regions
    )
    setattr(output, "class", c("vl_enrRank", "list"))
    
    # Save ----
    saveRDS(output, tmp)
    
  } else
    output <- readRDS(tmp)
   
  # Return ----
  return(output)
}