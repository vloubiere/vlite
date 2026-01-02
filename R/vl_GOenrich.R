#' Gene Ontology (GO) Enrichment Analysis
#'
#' This function computes Gene Ontology (GO) enrichment for a group of genes or clusters of genes.
#' It calculates enrichment statistics (log2 odds ratio and adjusted p-values) for selected GO annotations
#' (Biological Process, Cellular Component, and Molecular Function) and optionally plots the results.
#'
#' @param geneIDs A vector (for a single group) or a named list (for multiple clusters) of gene IDs to analyze.
#' @param geneUniverse.IDs A vector of gene IDs defining the universe of genes for the analysis. If `NULL`,
#' all genes are used as the universe.
#' @param species Character. The species for the analysis. Choose between `Dm` (Drosophila melanogaster)
#' and `Mm` (Mus musculus).
#' @param keyType The keyType to be used. Defaults to 'FLYBASE' (when genome= 'Dm') or 'ENSEMBL' (when genome= 'Mm').
#' For Dm, 'SYMBOL' can be used for gene symbols.
#' @param select Character vector specifying which GO annotations to consider. Options are `BP`
#' (Biological Process), `CC` (Cellular Component), and `MF` (Molecular Function). Default= c("BP", "CC", "MF").
#' @param log2OR.pseudocount Numeric. A pseudocount added to the contingency table to avoid infinite values
#' in the log2 odds ratio calculation. Default= 0.5.
#' @param cleanup.cache Logical. If `TRUE`, clears cached intermediate results. Default= FALSE.
#'
#' @details
#' The function performs GO enrichment analysis by comparing the overlap of input gene sets with GO terms
#' against a background universe of genes. It uses Fisher's exact test to compute p-values and calculates
#' log2 odds ratios for enrichment. Adjusted p-values (FDR) are computed for multiple testing correction.
#'
#' If multiple clusters of genes are provided, the function computes enrichment for each cluster separately.
#' The results are returned as a data table with class `vl_enr` (for a single group) or `vl_enr_cl`
#' (for multiple clusters), which can be directly plotted using the appropriate plotting methods.
#'
#' @return A data table containing the enrichment results with the following columns:
#' - `cl`: Cluster name (if multiple clusters are analyzed).
#' - `variable`: GO term ID.
#' - `annotation`: GO annotation type (e.g., BP, CC, MF).
#' - `name`: GO term description.
#' - `set_hit`: Number of genes in the input set associated with the GO term.
#' - `set_total`: Total number of genes in the input set.
#' - `universe_hit`: Number of genes in the universe associated with the GO term.
#' - `universe_total`: Total number of genes in the universe.
#' - `log2OR`: Log2 odds ratio for enrichment.
#' - `padj`: Adjusted p-value (FDR).
#'
#' @examples
#' # Example using ORFtag hits from the Nat. Methods paper -------
#' # Download data
#' tmp <- tempfile(pattern = ".xlsx")
#' download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-024-02339-x/MediaObjects/41592_2024_2339_MOESM3_ESM.xlsx",
#'               destfile = tmp)
#'
#' # Activators
#' act <- readxl::read_xlsx(tmp, sheet = 1)
#' act <- act$gene_id[act$hit==TRUE]
#'
#' # Repressors
#' rep <- readxl::read_xlsx(tmp, sheet = 2)
#' rep <- rep$gene_id[rep$hit==TRUE]
#' 
#' # Background (genes with input integration)
#' bg <- readxl::read_xlsx(tmp, sheet = 1)
#' bg <- bg$gene_id[!is.na(bg$FDR)]
#'
#' # Compute activators enrichment
#' enr <- vl_GOenrich(geneIDs = act,
#'                    geneUniverse.IDs = bg,
#'                    species= "Mm")
#'
#' # Plot
#' vl_par(mai= c(.9, 2, .9, 1.3))
#' plot(obj= enr[annotation=="CC"],
#'      top.enrich = 10,
#'      order= "log2OR")
#'
#' # Compare activators and repressor Cellular Compartments
#' enr2 <- vl_GOenrich(geneIDs = list(Act= act,
#'                                    Rep= rep),
#'                     geneUniverse.IDs = bg,
#'                     species= "Mm",
#'                     select= "CC")
#'
#' # Plot
#' vl_par(mai= c(.9, 2, .9, 1.3))
#' plot(obj= enr2,
#'      top.enrich = 10,
#'      padj.cutoff= 0.001,
#'      order = "log2OR",
#'      cex= .5)
#'
#' # Other example using gene clusters from the epigenetic cancer nature paper -------
#' # Download EpiCancer clusters
#' tmp <- tempfile(fileext = ".xlsx")
#' download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07328-w/MediaObjects/41586_2024_7328_MOESM5_ESM.xlsx",
#'               destfile = tmp)
#'
#' # Import genes
#' genes <- readxl::read_xlsx(tmp, sheet = 1)
#' genes <- as.data.table(genes)
#'
#' # Gene clusters
#' cl <- genes[!cluster %in% c("NA", "Unaffected")]
#' cl <- split(cl$FBgn, cl[, .(cluster, V2= ifelse(PcG_bound, "Bound", "Unbound"))])
#'
#' # Compare clusters
#' enr <- vl_GOenrich(geneIDs = cl,
#'                    geneUniverse.IDs = genes$FBgn,
#'                    species= "Dm")
#' enr[, cl:= factor(cl,
#'                   c("Reversible.Bound",
#'                     "Irreversible.Bound",
#'                     "Transient-specific.Bound",
#'                     "Down 1.Bound",
#'                     "Down 2.Bound",
#'                     "Down 3.Bound",
#'                     "Reversible.Unbound",
#'                     "Irreversible.Unbound",
#'                     "Transient-specific.Unbound",
#'                     "Down 1.Unbound",
#'                     "Down 2.Unbound",
#'                     "Down 3.Unbound"))]
#'
#' # Plot
#' vl_par(mai= c(.9, 2, .9, 1.3))
#' plot(obj= enr,
#'      top.enrich = 5,
#'      padj.cutoff= 0.05,
#'      order = "log2OR",
#'      cex= .5)
#'
#' @export
vl_GOenrich <- function(geneIDs,
                        geneUniverse.IDs,
                        species,
                        keyType,
                        select= c("BP", "CC", "MF"),
                        log2OR.pseudocount= 0.5,
                        cleanup.cache= FALSE)
{
  # Checks ----
  if(!species %in% c("Dm", "Mm"))
    stop("Species should be one of 'Dm' or 'Mm'")
  if(!is.character(geneUniverse.IDs))
    stop("geneUniverse.IDs should be a character vector of gene IDs.")
  # Format and make unique ----
  if(is.character(geneIDs))
    geneIDs <- list(set= geneIDs)
  if(is.null(names(geneIDs)))
    names(geneIDs) <- seq_along(geneIDs)
  if(anyDuplicated(names(geneIDs)))
    stop("names(geneIDs) are not unique.")
  geneIDs <- lapply(geneIDs, unique)
  geneUniverse.IDs <- unique(geneUniverse.IDs)
  if(any(!unlist(geneIDs) %in% geneUniverse.IDs))
    stop("Some geneIDs were absent from the geneUniverse.IDs.")
  
  # Genome ----
  db <- switch(species,
               "Dm"= org.Dm.eg.db::org.Dm.eg.db,
               "Mm"= org.Mm.eg.db::org.Mm.eg.db)
  if(missing(keyType))
    keyType <- switch(species,
                      "Dm"= "FLYBASE",
                      "Mm"= "ENSEMBL")
  
  # Temp directory for caching GOs ----
  cache.file <- vl_cache_file(list(geneUniverse.IDs, species))
  
  # Extract GOs associated to universe ----
  if(cleanup.cache | !file.exists(cache.file)) {
    set <- AnnotationDbi::select(x= db,
                                 keys = geneUniverse.IDs,
                                 keytype= keyType,
                                 columns= "GOALL")
    set <- data.table::as.data.table(set)
    # Simplify
    set$EVIDENCEALL <- NULL
    setnames(set,
             c("gene_id", "GO", "annotation"))
    set <- na.omit(unique(set))
    # Add description
    terms <- AnnotationDbi::select(
      GO.db::GO.db,
      keys = as.character(unique(set$GO)),
      keytype= "GOID",
      columns= "TERM"
    )
    terms <- as.data.table(terms)
    set[terms, name:= i.TERM, on= "GO==GOID"]
    # Count
    set[, universe_hit:= uniqueN(gene_id), .(GO, annotation)]
    set[, universe_total:= uniqueN(gene_id), annotation]
    # Save
    saveRDS(set, cache.file)
  } else
    set <- readRDS(cache.file)
  
  # Select annotations ----
  set <- set[annotation %in% select]
  
  # Compute overlaps ----
  enr <- lapply(geneIDs, function(x) data.table(gene_id= unique(x)))
  enr <- rbindlist(enr, idcol = "cl")
  enr <- merge(enr, set, by= "gene_id")
  enr[, set_hit:= uniqueN(gene_id), .(cl, GO, annotation)]
  enr[, set_total:= uniqueN(gene_id), .(cl, annotation)]
  
  # Warning genes with no annotations ----
  message("Genes mapped:")
  message(paste("Set:", uniqueN(enr$gene_id), "/", uniqueN(unlist(geneIDs))))
  message(paste("Universe:", uniqueN(set$gene_id), "/", uniqueN(geneUniverse.IDs)))
  
  # Compute enrichments per cluster ----
  enr[, c("OR", "pval"):= {
    # Contingency matrix
    mat <- c(
      set_hit,
      set_total - set_hit,
      universe_hit - set_hit,
      (universe_total - set_total) - (universe_hit - set_hit)
    )
    mat <- matrix(mat, byrow= T, ncol= 2)
    # p.value
    .f <- fisher.test(mat, alternative = "greater")
    # log2OR (pseudocount avoid Inf)
    if(any(mat==0)) {
      mat <- mat+log2OR.pseudocount
      .f$estimate <- (mat[1,1] * mat[2,2]) / (mat[2,1] * mat[1,2])
    }
    .(.f$estimate, .f$p.value)
  }, .(set_hit, set_total, universe_hit, universe_total)]
  
  # Make unique object ----
  uniq.cols <- setdiff(names(enr), "gene_id")
  res <- enr[, .(associated_genes= paste0(sort(unique(gene_id)), collapse = ",")), uniq.cols]
  
  # Compute log2OR and padj ----
  res[, log2OR:= log2(OR)]
  res[, padj:= p.adjust(pval, method = "fdr"), .(cl, annotation)]
  res$OR <- res$pval <- NULL
  
  # Order ----
  res[, cl:= factor(cl, names(geneIDs))]
  setorderv(res, c("cl", "padj"))
  setcolorder(res, c("cl", "GO", "annotation", "name", "log2OR", "padj"))
  
  # Define class (for plotting methods) ----
  if(length(geneIDs)>1) {
    setattr(res, "class", c("vl_enr_cl", "data.table", "data.frame"))
  } else {
    setattr(res, "class", c("vl_enr", "data.table", "data.frame"))
  }
  
  # Order and return enrichment table ----
  return(res)
}
