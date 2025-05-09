% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_vl_enr.R
\name{plot.vl_enr}
\alias{plot.vl_enr}
\title{Plot Enrichment Results for vl_enr Objects}
\usage{
\method{plot}{vl_enr}(
  obj,
  padj.cutoff = 0.05,
  top.enrich = Inf,
  min.counts = 3L,
  order = "log2OR",
  xlab = "Odd Ratio (log2)",
  breaks = NULL,
  col = c("blue", "red")
)
}
\arguments{
\item{obj}{An object of class `vl_enr`. This should be a data table containing enrichment results,
including columns for `padj` (adjusted p-values), `log2OR` (log2 odds ratio), and `set_hit` (counts).}

\item{padj.cutoff}{Numeric. The adjusted p-value cutoff to filter enrichments before plotting.
Default is `0.05`.}

\item{top.enrich}{Integer. The maximum number of top enrichments to plot, based on the `order` parameter.
Default is `Inf` (plot all enrichments that pass the cutoffs).}

\item{min.counts}{Integer. The minimum number of counts (`set_hit`) required to include an enrichment
in the plot. Default is `3L`.}

\item{order}{Character. The metric used to order enrichments before selecting the top enrichments.
Possible values are `padj` (adjusted p-value) or `log2OR` (log2 odds ratio). Default is `log2OR`.}

\item{xlab}{Character. The label for the x-axis of the barplot. Default is `Odd Ratio (log2)`.}

\item{breaks}{Numeric vector. The color breaks to use for the color scale. Defaults to the range of
filtered `padj` values (-log10 transformed).}

\item{col}{Character vector. The color scale to use for the barplot. Default is `c("blue", "red")`.}
}
\value{
Invisibly returns the filtered and plotted data table.
}
\description{
This function generates a barplot to visualize enrichment results from an object of class `vl_enr`.
The bars represent the log2 odds ratio (log2OR) of enrichment, and the color scale corresponds to
the adjusted p-values (-log10 transformed).
}
\details{
The function filters the input data based on the specified `padj.cutoff` and `min.counts` thresholds.
It then orders the data based on the `order` parameter and selects the top enrichments to plot.
The barplot displays the log2 odds ratio (log2OR) for each enrichment, with colors representing
the adjusted p-values (-log10 transformed). A heat key is also displayed to indicate the color scale.
}
\examples{
# Example usage:
For an example using GO enrichment, see ?vl_GOenrich()
For an example using motif enrichment, see ?vl_motifEnrich()

# Assuming `enrichment_results` is an object of class `vl_enr`:
plot.vl_enr(enrichment_results, padj.cutoff = 0.01, top.enrich = 10, min.counts = 5)

}
