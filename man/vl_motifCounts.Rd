% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vl_motifCounts.R
\name{vl_motifCounts}
\alias{vl_motifCounts}
\title{Count Motif Occurrences in Genomic Regions}
\usage{
vl_motifCounts(
  sequences,
  bed,
  pwm_log_odds,
  genome,
  bg = "genome",
  p.cutoff = 0.00005,
  cleanup.cache = FALSE
)
}
\arguments{
\item{sequences}{A named character vector of sequences to analyze. This argument takes precedence over the bed argument.}

\item{bed}{Genomic ranges in a format compatible with ?importBed, from which genomic sequences will be retrieved when sequences is set to NULL.}

\item{pwm_log_odds}{A PWMatrixList (in log2 odds ratio format) containing motifs to count.}

\item{genome}{Genome to use as background when bg = "genome" and/or to retrieve sequences (when bed is specified).}

\item{bg}{Background model for motif detection. Options are "genome", "subject" (inferred from input sequences) or "even" (0.25, 0.25, 0.25, 0.25). Default= "genome".}

\item{p.cutoff}{p-value cutoff for motif detection. Default= 5e-5.}

\item{cleanup.cache}{Logical. If set to TRUE, clears cached intermediate results. Default= FALSE.}
}
\value{
A matrix of motif counts.
}
\description{
Counts motif occurrences in a set of sequences or genomic regions using a PWMatrixList.
}
\examples{
# Download Dev enhancer from pe-STARR-Seq paper
tmp <- tempfile(pattern = ".xlsx")
download.file(url = "https://static-content.springer.com/esm/art\%3A10.1038\%2Fs41467-024-52921-2/MediaObjects/41467_2024_52921_MOESM4_ESM.xlsx",
              destfile = tmp)

# Retrieve enhancers
dev <- readxl::read_xlsx(tmp, sheet = 2, skip = 1)
dev <- as.data.table(dev)
dev <- dev[group=="dev" & detail \%in\% c("medium", "strong")]

# Load motifs and select JASPAR
load("/groups/stark/vloubiere/motifs_db/vl_Dmel_motifs_DB_full.RData")
sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds]

# Compute counts from regions
mot1 <- vl_motifCounts(bed= dev,
                       genome= "dm3",
                       pwm_log_odds= sel)

# Compute counts from sequences
mot2 <- vl_motifCounts(sequences = dev$enhancer_sequence,
                       genome= "dm3",
                       pwm_log_odds= sel)

# Both approaches should be identical
identical(mot1, mot2)

}
