% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_vl_enr_clusters.R
\name{plot.vl_enr_cl}
\alias{plot.vl_enr_cl}
\title{Plot Clustered Enrichment Results for vl_enr_cl Objects}
\usage{
\method{plot}{vl_enr_cl}(
  obj,
  padj.cutoff = 0.05,
  top.enrich = Inf,
  min.counts = 3L,
  order = "log2OR",
  color.breaks = NULL,
  size.legend.breaks = NULL,
  cex = 1,
  col = c("blue", "red"),
  main = NA,
  plot.empty.clusters = T
)
}
\arguments{
\item{obj}{An object of class `vl_enr_cl`. This should be a data table containing clustered enrichment results,
including columns for `padj` (adjusted p-values), `log2OR` (log2 odds ratio), `set_hit` (counts),
`name` (feature names), and `cl` (clusters).}

\item{padj.cutoff}{Numeric. The adjusted p-value cutoff to filter enrichments before plotting.
Default is `0.05`.}

\item{top.enrich}{Integer. The maximum number of top enrichments to plot per cluster, based on the `order` parameter.
Default is `Inf` (plot all enrichments that pass the cutoffs).}

\item{min.counts}{Integer. The minimum number of counts (`set_hit`) required to include an enrichment
in the plot. Default is `3L`.}

\item{order}{Character. The metric used to order enrichments before selecting the top enrichments.
Possible values are `padj` (adjusted p-value) or `log2OR` (log2 odds ratio). Default is `log2OR`.}

\item{color.breaks}{Numeric vector. The color breaks to use for the color scale. Defaults to the range of
filtered `padj` values (-log10 transformed).}

\item{size.legend.breaks}{Numeric vector. Breakpoints for the size legend in the balloon plot. Default is `NULL`.}

\item{cex}{Numeric. Scaling factor for the plot. Default is `1`.}

\item{col}{Character vector. The color scale to use for the balloon plot. Default is `c("blue", "red")`.}

\item{main}{Character. The main title for the plot. Default is `NA` (no title).}

\item{plot.empty.clusters}{Logical. Whether to include clusters with no enrichments in the plot.
Default is `TRUE`.}
}
\value{
Invisibly returns the filtered and plotted data table.
}
\description{
This function generates a balloon plot to visualize clustered enrichment results from an object of class `vl_enr_cl`.
The balloons represent the log2 odds ratio (log2OR) of enrichment, with their size corresponding to the odds ratio
and their color representing the adjusted p-values (-log10 transformed).
}
\details{
The function filters the input data based on the specified `padj.cutoff`, `min.counts`, and `log2OR > 0` thresholds.
It then orders the data based on the `order` parameter and selects the top enrichments per cluster.
The balloon plot displays the log2 odds ratio (log2OR) as the size of the balloons, with colors representing
the adjusted p-values (-log10 transformed). Clusters with no enrichments can be optionally excluded from the plot.
}
\examples{
# Example usage:
For an example using GO enrichment, see ?vl_GOenrich()
For an example using motif enrichment, see ?vl_motifEnrich()

# Assuming `clustered_enrichment_results` is an object of class `vl_enr_cl`:
plot.vl_enr_cl(clustered_enrichment_results, padj.cutoff = 0.01, top.enrich = 5, min.counts = 5)

}
