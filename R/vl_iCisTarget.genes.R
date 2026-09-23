# Motif enrichment analysis using rCisTarget
#' 
#' A wrapper around the rCisTarget package that identify motif enriched in Dmel genes.
#'
#' @param geneList A named list of Dmel gene symbols.
#' @param highlightTFs If a character list of TF names is provided, the column TFinDB in the otuput table will indicate whether
#' any of those TFs are included within the 'high-confidence' annotation (two asterisks, **) or 'low-confidence' annotation
#' (one asterisk, *) of the motif. The vector can be named to indicate which TF to highlight for each gene-set. 
#' Otherwise, all TFs will be used for all geneSets.
#' @param NES.cutoff NES threshold to identify significant motifs. The NES is calculated -for each motif- based on the AUC distribution of all the motifs for the gene-set [(x-mean)/sd]. Default= 3.0.
#' @param maxRank The maximum rank to take into account for the gene enrichment recovery curve. Default= 5000.
#' @param aucMaxRank Threshold to calculate the AUC. In a simplified way, the AUC value represents the fraction of genes
#' (within the top X genes in the ranking) that are included in the signature. Default= 0.05 (5% of the total number of
#' genes in the rankings). Common values: 1-10%.
#' @param plot Should the result be plotted?
#' @param cluster.rows Should the rows of the resulting heatmap be clustered? Default= FALSE
#' @param cluster.cols Should the columns of the resulting heatmap be clustered? Default= FALSE
#' @param drop.empty Should sub-lists with no enrichment be dropped? Default= FALSE.
#' @param legend.title Default= "NES".
#' @param show.grid Should the grid be plotted? Default= TRUE.
#' @param show.numbers Should the NES values be plotted? Default= TRUE.
#' @param motifRankings.feather Path to the gene motifRankings .feather file. Default=
#' "/zssd/scratch/vincent.loubiere/motifs_db/iCisTarget/mc_v10_clust/gene_based/dm6_v10_clust.genes_vs_motifs.rankings.feather".
#' @param motifAnnot.tbl Path to the gene motifAnnot .tbl file. Default=
#' "/zssd/scratch/vincent.loubiere/motifs_db/iCisTarget/motif2tf/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"
#' @param motifAnnot.feather 
#' @param nCores 
#' @param ... Extra arguments passed to the vl_heatmap function. 
#'
#' @returns
#' @export
#'
#' @examples
vl_iCisTarget.genes <- function(
    geneList,
    highlightTFs= NULL,
    NES.cutoff= 3.0,
    maxRank= 5000,
    aucMaxRank= 0.05,
    plot= T,
    cluster.rows= F,
    cluster.cols= F,
    drop.empty= FALSE,
    legend.title= "NES",
    show.grid= T,
    show.numbers= T,
    motifRankings.feather=
      "/zssd/scratch/vincent.loubiere/motifs_db/iCisTarget/mc_v10_clust/gene_based/dm6_v10_clust.genes_vs_motifs.rankings.feather",
    motifAnnot.tbl=
      "/zssd/scratch/vincent.loubiere/motifs_db/iCisTarget/motif2tf/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl",
    nCores= 4,
    cleanup.cache = FALSE,
    ...
)
{
  # Checks
  stopifnot(!is.null(names(geneList)))
  
  # Check if already computed
  output.file <- vl_cache_file(
    input.list = 
      list(
        geneList,
        highlightTFs,
        maxRank,
        aucMaxRank,
        motifRankings.feather,
        motifAnnot.tbl,
        nCores
      )
  )
  
  # If not compute yet
  if(cleanup.cache || !file.exists(output.file)) {
    
    # Import databases
    motifRankings <- RcisTarget::importRankings(motifRankings.feather)
    motifAnnot <- if(is.null(motifAnnot.tbl)) NULL else RcisTarget::importAnnotations(motifAnnot.tbl)
    
    # Compute enrichment ----
    res <- RcisTarget::cisTarget(
      geneSets = geneList,
      motifRankings = motifRankings,
      motifAnnot = motifAnnot,
      motifAnnot_highConfCat = c("directAnnotation", "inferredBy_Orthology"),
      motifAnnot_lowConfCat = c("inferredBy_MotifSimilarity",
                                "inferredBy_MotifSimilarity_n_Orthology"),
      highlightTFs = NULL,
      nesThreshold = 3.0, # The cutoff will be done later (to avoid many re-runs)
      aucMaxRank = aucMaxRank * motifRankings@nColsInDB,
      geneErnMethod = "aprox",
      geneErnMaxRank = maxRank,
      nCores = nCores,
      verbose = TRUE
    )
    
    # Keep geneSet order ----
    res[, geneSet:= factor(geneSet, names(geneList))]
    
    # TF simplified names ----
    res[, TF_name:= gsub(" \\(.*$", "", TF_highConf)]
    res[, TF_name:= paste0(sort(unique(unlist(tstrsplit(TF_name, "; ")))), collapse = ";"), TF_name]
    
    # Save ----
    saveRDS(res, output.file)
  } else
    res <- readRDS(output.file)
  
  # NES cutoff ----
  res <- res[NES>=NES.cutoff]
  
  # Adjust TF order after NES cutoff ----
  setorderv(res, c("geneSet", "NES"), c(1, -1))
  res[, TF_name:= factor(TF_name, unique(TF_name))]
  
  # Plot ----
  if(plot) {
    # Dcast max NES
    mat <- dcast(res, TF_name~geneSet, value.var = "NES", fun.aggregate = function(x) max(c(0, x)), drop = drop.empty)
    mat <- as.matrix(mat, 1)
    # Extract numbers
    if(isTRUE(show.numbers)) {
      show.numbers <- round(mat, 1)
      show.numbers[show.numbers==0] <- NA
    }
    # Plot heatmap
    vl_heatmap(
      mat,
      cluster.rows = cluster.rows,
      cluster.cols = cluster.cols,
      show.grid = show.grid,
      show.numbers = show.numbers,
      legend.title = legend.title,
      ...
    )
  }
  
  # Return object ----
  return(res)
}