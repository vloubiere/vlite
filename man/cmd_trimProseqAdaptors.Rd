% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmd_trimProseqAdaptors.R
\name{cmd_trimProseqAdaptors}
\alias{cmd_trimProseqAdaptors}
\title{Generate Cutadapt Commands for Trimming PRO-Seq Adapters}
\usage{
cmd_trimProseqAdaptors(fq1, fq.output.folder = "db/fq/")
}
\arguments{
\item{fq1}{A character vector of .fq (or .fq.gz) file paths.}

\item{fq.output.folder}{Directory for trimmed FASTQ files. Default: `"db/fq/"`.}
}
\value{
A `data.table` with:
- `file.type`: Output file label.
- `path`: Path to the trimmed FASTQ file.
- `cmd`: Shell command to run Cutadapt.
}
\description{
Creates shell commands to trim PRO-Seq adapters from gzipped FASTQ files using Cutadapt.
}
\examples{
cmd <- cmd_trimProseqAdaptors(fq1 = "sample_R1.fq.gz", fq.output.folder = "db/fq/")
vl_submit(cmd, execute= FALSE)

}
