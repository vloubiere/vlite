% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contribSeqLogo.R
\name{contribSeqLogo}
\alias{contribSeqLogo}
\title{Plot contribution scores matrix}
\usage{
contribSeqLogo(contrib, row.idx = 1, start, end, highlight.pos)
}
\arguments{
\item{contrib}{A contrib object imported using ?importContrib.}

\item{row.idx}{The line index within the contrib object that should be plotted. Default= 1.}

\item{start}{Integer specifying at which the start position of the sequence should be clipped.}

\item{end}{Integer specifying at which the end position of the sequence should be clipped.}

\item{higlight}{Index of the positions that should be higlighted.}
}
\description{
Plot contribution scores matrix
}
