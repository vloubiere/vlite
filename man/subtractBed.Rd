% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bedTools.R
\name{subtractBed}
\alias{subtractBed}
\title{Subtract Genomic Regions from Reference Intervals}
\usage{
subtractBed(a, b, minoverlap = 1L, ignore.strand = TRUE)
}
\arguments{
\item{a}{Query regions in any format compatible with ?importBed.}

\item{b}{Target regions in any format compatible with ?importBed.}

\item{minoverlap}{A single integer specifying the minimum overlap with regions in for them to
be subtracted. Default= 1L.}

\item{ignore.strand}{If set to FALSE, only subtracts features that are on the same strand.
If set to TRUE (default), subtracts overlapping feature(s) regardless of their strand(s).}
}
\value{
A data.table containing the remaining portions of a, after subtracting the regions
defined in b.
}
\description{
A wrapper around ?GenomicRanges::setdiff that subtracts the regions in b to the genomic ranges in a.
}
\examples{
# Create example regions
a <- importBed(c("chr1:200-300:-", "chr1:100-500:-"))
b <- importBed(c("chr1:200-300:+", "chr1:400-450:-", "chr1:425-475:-"))

# Basic example
subtractBed(a, b)

# Only subtract regions with similar strand
subtractBed(a, b, ignore.strand= FALSE)

}
