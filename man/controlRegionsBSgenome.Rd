% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/controlRegionsBSgenome.R
\name{controlRegionsBSgenome}
\alias{controlRegionsBSgenome}
\title{Generate control regions}
\usage{
controlRegionsBSgenome(bed, genome, no.overlap = FALSE)
}
\arguments{
\item{bed}{Bed file used to produce similar control}

\item{genome}{BSgneome object to use. ex: "dm3", "dm6"}

\item{no.overlap}{If set to TRUE, avoids overlap between control sequences and the original bed file. Default= FALSE.}
}
\value{
data.table containing control regions
}
\description{
Generate regions with similar dispersion, seqnames and widths distribution than the provided bed file
}
\examples{
controlRegionsBSgenome(vl_SUHW_top_peaks, "dm3")

}
