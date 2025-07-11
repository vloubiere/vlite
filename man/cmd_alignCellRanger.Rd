% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmd_alignCellRanger.R
\name{cmd_alignCellRanger}
\alias{cmd_alignCellRanger}
\title{Generate Cell Ranger Alignment Command}
\usage{
cmd_alignCellRanger(
  fq.prefix,
  index,
  output.prefix,
  output.folder = "db/scRNASeq_10X/",
  cores = 10,
  mem = 64
)
}
\arguments{
\item{fq.prefix}{Prefix shared between input fastq files: '/path/to/fastq/shared_prefix'.}

\item{index}{Path to the Cell Ranger reference transcriptome. e.g.: '/groups/stark/vloubiere/genomes/Drosophila_melanogaster/cellranger_dm6_TG/index/'}

\item{output.prefix}{Prefix for output files.}

\item{output.folder}{Directory where output files will be saved. Default= 'db/scRNASeq_10X/'.}

\item{cores}{Number of CPU cores to use. Default= 10.}

\item{mem}{Memory in GB to allocate. Default= 64.}
}
\description{
Constructs a shell command to run 10x Genomics Cell Ranger for single-cell RNA-seq alignment and quantification, with resource specification and output management.
}
\examples{
cmd_alignCellRanger(
  fq.prefix = "/data/fastq/sample1",
  index = "/data/genomes/cellranger_index",
  output.prefix = "sample1",
  output.folder = "db/scRNASeq_10X/"
)

}
