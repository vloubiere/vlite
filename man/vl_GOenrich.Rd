% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vl_GOenrich.R
\name{vl_GOenrich}
\alias{vl_GOenrich}
\title{Gene Ontology (GO) Enrichment Analysis}
\usage{
vl_GOenrich(
  geneIDs,
  geneUniverse.IDs,
  species,
  select = c("BP", "CC", "MF"),
  log2OR.pseudocount = 1,
  cleanup.cache = FALSE
)
}
\arguments{
\item{geneIDs}{A vector (for a single group) or a named list (for multiple clusters) of gene IDs to analyze.}

\item{geneUniverse.IDs}{A vector of gene IDs defining the universe of genes for the analysis. If `NULL`,
all genes are used as the universe.}

\item{species}{Character. The species for the analysis. Choose between `Dm` (Drosophila melanogaster)
and `Mm` (Mus musculus).}

\item{select}{Character vector specifying which GO annotations to consider. Options are `BP`
(Biological Process), `CC` (Cellular Component), and `MF` (Molecular Function). Default= c("BP", "CC", "MF").}

\item{log2OR.pseudocount}{Numeric. A pseudocount added to the contingency table to avoid infinite values
in the log2 odds ratio calculation. Default= 1L.}

\item{cleanup.cache}{Logical. If `TRUE`, clears cached intermediate results. Default= FALSE.}
}
\value{
A data table containing the enrichment results with the following columns:
- `cl`: Cluster name (if multiple clusters are analyzed).
- `variable`: GO term ID.
- `annotation`: GO annotation type (e.g., BP, CC, MF).
- `name`: GO term description.
- `set_hit`: Number of genes in the input set associated with the GO term.
- `set_total`: Total number of genes in the input set.
- `universe_hit`: Number of genes in the universe associated with the GO term.
- `universe_total`: Total number of genes in the universe.
- `log2OR`: Log2 odds ratio for enrichment.
- `padj`: Adjusted p-value (FDR).
}
\description{
This function computes Gene Ontology (GO) enrichment for a group of genes or clusters of genes.
It calculates enrichment statistics (log2 odds ratio and adjusted p-values) for selected GO annotations
(Biological Process, Cellular Component, and Molecular Function) and optionally plots the results.
}
\details{
The function performs GO enrichment analysis by comparing the overlap of input gene sets with GO terms
against a background universe of genes. It uses Fisher's exact test to compute p-values and calculates
log2 odds ratios for enrichment. Adjusted p-values (FDR) are computed for multiple testing correction.

If multiple clusters of genes are provided, the function computes enrichment for each cluster separately.
The results are returned as a data table with class `vl_enr` (for a single group) or `vl_enr_cl`
(for multiple clusters), which can be directly plotted using the appropriate plotting methods.
}
\examples{
# Example using ORFtag hits from the Nat. Methods paper -------
# Download data
tmp <- tempfile(pattern = ".xlsx")
download.file(url = "https://static-content.springer.com/esm/art\%3A10.1038\%2Fs41592-024-02339-x/MediaObjects/41592_2024_2339_MOESM3_ESM.xlsx",
              destfile = tmp)

# Activators
act <- readxl::read_xlsx(tmp, sheet = 1)
act <- act$gene_id[act$hit==TRUE]

# Background (genes with input integration)
bg <- readxl::read_xlsx(tmp, sheet = 1)
bg <- bg$gene_id[!is.na(bg$FDR)]

# Compute activators enrichment
enr <- vl_GOenrich(geneIDs = act,
                   geneUniverse.IDs = bg,
                   species= "Mm")

# Repressors
rep <- readxl::read_xlsx(tmp, sheet = 2)
rep <- rep$gene_id[rep$hit==TRUE]

# Plot
vl_par(mai= c(.9, 2, .9, 1.3))
plot(obj= enr[annotation=="CC"],
     top.enrich = 10,
     order= "log2OR")

# Compare activators and repressor Cellular Compartments
enr2 <- vl_GOenrich(geneIDs = list(Act= act,
                                   Rep= rep),
                    geneUniverse.IDs = bg,
                    species= "Mm",
                    select= "CC")

# Plot
vl_par(mai= c(.9, 2, .9, 1.3))
plot(obj= enr2,
     top.enrich = 10,
     padj.cutoff= 0.001,
     order = "log2OR",
     cex= .5)

# Other example using gene clusters from the epigenetic cancer nature paper -------
# Download EpiCancer clusters
tmp <- tempfile(fileext = ".xlsx")
download.file(url = "https://static-content.springer.com/esm/art\%3A10.1038\%2Fs41586-024-07328-w/MediaObjects/41586_2024_7328_MOESM5_ESM.xlsx",
              destfile = tmp)

# Import genes
genes <- readxl::read_xlsx(tmp, sheet = 1)
genes <- as.data.table(genes)

# Gene clusters
cl <- genes[!cluster \%in\% c("NA", "Unaffected")]
cl <- split(cl$FBgn, cl[, .(cluster, V2= ifelse(PcG_bound, "Bound", "Unbound"))])

# Compare activators and repressor Cellular Compartments
enr <- vl_GOenrich(geneIDs = cl,
                   geneUniverse.IDs = genes$FBgn,
                   species= "Dm")
enr[, cl:= factor(cl,
                  c("Reversible.Bound",
                    "Irreversible.Bound",
                    "Transient-specific.Bound",
                    "Down 1.Bound",
                    "Down 2.Bound",
                    "Down 3.Bound",
                    "Reversible.Unbound",
                    "Irreversible.Unbound",
                    "Transient-specific.Unbound",
                    "Down 1.Unbound",
                    "Down 2.Unbound",
                    "Down 3.Unbound"))]

# Plot
vl_par(mai= c(.9, 2, .9, 1.3))
plot(obj= enr,
     top.enrich = 5,
     padj.cutoff= 0.05,
     order = "log2OR",
     cex= .5)

}
