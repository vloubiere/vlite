% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importBed.R
\name{importBed}
\alias{importBed}
\title{Import Genomic Ranges in BED Format}
\usage{
importBed(bed, col.names = NULL)
}
\arguments{
\item{bed}{Input genomic ranges in one of the following formats:
\itemize{
  \item Character vector of genomic coordinates ("chr:start-end[:strand]").
  \item File path(s) to .bed, .narrowPeak, or .broadPeak files. If several files are provided, they will be catenated.
  \item GRanges object.
  \item data.frame or data.table with required columns.
}}

\item{col.names}{Can be used to overwrite default column names.}
}
\description{
Imports genomic coordinates from various input formats into a standardized gr data.table format.
Supported input types: GRanges objects, character strings ("chr2L:30000-31000:-"), .bed, .narrowPeak
and .broadPeak file paths (specifications can be found at https://genome.ucsc.edu/FAQ/FAQformat.html).
}
\details{
.bed, .narrowPeak and .broadPeak file specifications can be found at https://genome.ucsc.edu/FAQ/FAQformat.html.
Character coordinates should follow the syntax "chr:start-end[:strand]" (e.g. "chr2L:1000-2000:+", "chr2L:1000-2000"...).
GRanges objects are directly converted to data.table and metadata columns are preserved.
Required columns for data.frames/data.tables are 'seqnames', 'start', 'end' ('strand' is optional).
}
\examples{
# Import from character coordinates
importBed(c("chr2L:1000-2000:+", "chr3L:1000-2000:-", "chr3L:4000-6000:*"))[]

# Import from GRanges
gr <- GRanges("chr2L", IRanges(c(1000, 5000), c(2000, 6000)))
importBed(gr)[]

# Import single BED file
bed_file <- tempfile(fileext = ".bed")
rtracklayer::export(gr, bed_file)
importBed(bed_file)[]

# Import narrowpeak peak file
# Download example file
file1 <- tempfile(fileext = ".narrowPeak")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE222193&format=file&file#' =GSE222193\%5FCUTandRUN\%5FH3K27Ac\%5FControl\%5Fno\%5Fph\%2DKD\%5Fconfident\%5Fpeaks\%2EnarrowPeak\%2Egz",
              destfile = file1)
# Import
importBed(bed = file1)[]

# Import and merge multiple peak files               
# Download second example file
file2 <- tempfile(fileext = ".narrowPeak")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE222193&format=file&file#' =GSE222193\%5FChIP\%5FPH\%5FControl\%5Fno\%5Fph\%2DKD\%5Fconfident\%5Fpeaks\%2EnarrowPeak\%2Egz",
              destfile = file2)
# Import both at once (combined)
importBed(bed = c(file1, file2))[]

}
