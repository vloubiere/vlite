#' Title
#'
#' @param regions A list of regions to be overlapped with db peaks to compute enrichment.
#' @param max.rank Default= 5000
#' @param N.permutations Default= NA (no pvalues).
#' @param decreasing Default= FALSE.
#' @param seed Default= 1.
#' @param cleanup.cache Default= FALSE.
#' @param ranks.file The ranks to be used. Default to max cluster buster scores and number of cluster motifs,
#' saved at "/zssd/scratch/vincent.loubiere/projects/iCis_lite/Rdata/ranks/v2_core/ranks2.rds".
#'
#' @returns
#' @export
#'
#' @examples
vl_iCisEnrich <- function(
    regions,
    max.rank= 5000,
    N.permutations= NA,
    decreasing= F, 
    seed= 1,
    cleanup.cache = FALSE,
    what= "enr",
    pwm.db.file= "/zssd/scratch/vincent.loubiere/motifs_db/JASPAR_CORE/20260812_JASPAR2026_CORE_NON_REDUNDANT.rds",
    ranks.file= "/zssd/scratch/vincent.loubiere/projects/iCis_lite/Rdata/ranks/v2_core/ranks2.rds",
    tmp.dir= "/zssd/scratch/vincent.loubiere/projects/iCis_lite/db/tmp/"
) {
  # Checks
  if(is.data.table(regions) | !is.list(regions))
    regions <- list(set= regions)
  stopifnot(!is.null(names(regions)))
  what <- match.arg(arg = what, c("enr", "recov", "targets"), several.ok = T)
  
  # Check if exists
  tmp.file <- vlite::vl_cache_file(
    input.list = list(
      regions,
      max.rank,
      N.permutations,
      decreasing,
      seed,
      cleanup.cache,
      pwm.db.file,
      ranks.file,
      tmp.dir
    ),
    tmp.dir = tmp.dir
  )
  if(!file.exists(tmp.file) | cleanup.cache) {
    # Import regions
    regions <- lapply(regions, importBed)
    
    # Import ranks and peaks
    ranks <- readRDS(ranks.file)
    peaks <- ranks[, 1:3]
    ranks <- ranks[, -c(1:3)]
    
    # Compute enrichment
    enr <- lapply(regions, function(x) {
      vlite::vl_motifEnrichRanks(
        labels= covBed(peaks, x)>0,
        ranks.matrix= ranks,
        max.rank= max.rank,
        N.permutations= N.permutations,
        decreasing= decreasing,
        seed= seed,
        cleanup.cache= cleanup.cache,
        tmp.dir= tmp.dir
      )
    }
    )
    names(enr) <- names(regions)
    
    # Aggregate
    res <- list(
      enr= rbindlist(lapply(enr, `[[`, "enr"), idcol = "group"),
      targets= rbindlist(lapply(enr, `[[`, "candidate_regions"), idcol = "group"),
      recov= lapply(enr, `[[`, "recov")
    )
    # Annotate
    res$targets[, c("seqnames", "start", "end"):= peaks[res$targets$row.idx]]
    pwm.meta <- readRDS(pwm.db.file)$metadata[, .(motif= name, cluster, cluster.name)]
    res$enr <- merge(
      res$enr,
      pwm.meta,
      by= "motif"
    )
    setattr(res$enr, "class", c("vl_iCisEnr", "data.table"))
    setattr(res$recov, "class", c("vl_iCisRecov", "list"))
    
    # Save
    saveRDS(res, tmp.file)
  } else
    res <- readRDS(tmp.file)
   
  # Return
  if(length(what)==1)
    return(res[[what]]) else
      return(res[what])
}