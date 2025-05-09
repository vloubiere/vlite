% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binBSgenome.R
\name{binBSgenome}
\alias{binBSgenome}
\title{Bin BS genome}
\usage{
binBSgenome(
  genome,
  bins.width = 50L,
  steps.width = bins.width,
  restrict.seqnames = NULL
)
}
\arguments{
\item{genome}{BSgenome object to use. Example "dm3"}

\item{bins.width}{bins width default to 50bp}

\item{steps.width}{steps width separating each bin. default set to bins.width}

\item{restrict.seqnames}{If specified, bins are restricted to provided seqnames}
}
\value{
data.table containing bin coordinates
}
\description{
Very fast using data.table
}
\examples{
vl_binBSgenome(genome= "dm3",
               bins.width= 50)

}
