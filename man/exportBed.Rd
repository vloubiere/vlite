% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exportBed.R
\name{exportBed}
\alias{exportBed}
\title{Export Genomic Ranges to BED Format Files}
\usage{
exportBed(bed, file)
}
\arguments{
\item{bed}{Input genomic ranges, in any format compatible with ?importBed().}

\item{file}{Output file path ending with .bed, .narrowPeak or .broadPeak extensions.}
}
\description{
Exports genomic ranges as .bed, .narrowPeak or .broadPeak files
File specifications can be found at https://genome.ucsc.edu/FAQ/FAQformat.html.
}
\details{
Upon saving, 1-based start coordinates will be converted to 0-based, following BED format specifications.
}
\examples{
# Export simple BED format:
bed <- data.table(seqnames = "chr2L",
                  start = 1000,
                  end = 2000,
                  strand = "+")
exportBed(bed, file = "test.bed")
exportBed("chr3R:1000-2000:+", file = "test.bed")

# Export narrowPeak format:
peak_data <- data.table(
  seqnames = "chr2L",
  start = 1000,
  end = 2000,
  name = "peak1",
  score = 100,
  strand = "+",
  signalValue = 5.5,
  pValue = 0.001,
  qValue = 0.05,
  peak = 1500
)
exportBed(peak_data, file = "test.narrowPeak")

}
