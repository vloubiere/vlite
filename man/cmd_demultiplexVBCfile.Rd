% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmd_demultiplexVBCfile.R
\name{cmd_demultiplexVBCfile}
\alias{cmd_demultiplexVBCfile}
\title{Generate Demultiplexing Commands for .bam or .tar.gz Files from the VBC NGS Facility}
\usage{
cmd_demultiplexVBCfile(
  vbcFile,
  layout,
  i7 = "none",
  i5 = "none",
  umi = FALSE,
  i7.column = 14,
  i5.column = 12,
  output.prefix,
  fq.output.folder = "db/fq/",
  start.seq = NULL,
  trim.length,
  cores = 8,
  head
)
}
\arguments{
\item{vbcFile}{Path to the input BAM file or .tar.gz file containing the reads.}

\item{layout}{Sequencing layout, must be either "PAIRED" or "SINGLE".}

\item{i7}{i7 index sequence(s). Multiple indexes should be comma-separated.
If set to "none" (default), no i7 filtering is performed}

\item{i5}{i5 index sequence(s). Multiple indexes should be comma-separated.
If set to "none" (default), no i5 filtering is performed}

\item{i7.column}{Column number in the BAM file containing the i7 index. Default= 14L.}

\item{i5.column}{Column number in the BAM file containing the i5 index. Default= 12L.}

\item{output.prefix}{Output files prefix. If not provided,
constructed from the input file name and index sequences.}

\item{fq.output.folder}{Directory where output FASTQ files will be written. Default= "db/fq/".}

\item{start.seq}{For PRO-Seq data, the eBC DNA string that must be present at the start of the reads.
Only supported for BAM input. Default= NULL.}

\item{trim.length}{When start.seq is provided (PRO-Seq reads), number of nucleotides that should be cut
from the sequence (after start.seq has been trimmed) and appended to the read ID. Only supported for BAM input.}

\item{cores}{Number of CPU cores to use for samtools processing (when using BAM input). Default= 8L.}

\item{head}{Number of reads that should be processed (for testing purposes).}
}
\value{
A data.table with three columns:
\itemize{
  \item `file.type`: Labels for output files (e.g.: "fq1", "fq2").
  \item `path`: Full paths to the output files.
  \item `cmd`: Shell command(s) to generate the output files.
}
}
\description{
Generates shell commands to demultiplex BAM or tar.gz FASTQ files from the VBC NGS facility.
This function is a wrapper around the Perl scripts vbc_bam_demultiplexing.pl and vbc_tar_demultiplexing.pl.
It supports both single-end and paired-end sequencing data, with optional support for PRO-Seq-specific processing
(only supported with BAM input).
}
\details{
The function generates commands that:
1. Read a BAM file using samtools or a .tar.gz file.
2. Process the reads using vbc_bam_demultiplexing.pl (BAM) or vbc_tar_demultiplexing.pl (tar.gz).
3. Filter reads based on i7 and/or i5 indexes.
4. Optionally process PRO-Seq-specific requirements (BAM only).
5. Output gzipped FASTQ files.
}
\section{Output Files}{

For paired-end data:
- `<output.prefix>_1.fq.gz`
- `<output.prefix>_2.fq.gz`

For single-end data:
- `<output.prefix>.fq.gz`
}

\examples{
# Example with a .tar.gz file (One of Tomas screens)
cmd <- cmd_demultiplexVBCfile(
  vbcFile = "/groups/stark/projects/PE75_20250120_TP/2415402230_1_R18342_20250120.tar.gz",
  layout = "PAIRED",
  i7 = "CTATAC",
  i5 = "none",
  output.prefix = "ORFeome_fastq_demult_CTATAC",
  fq.output.folder = "/groups/stark/vloubiere/packages/tests/",
  head = 40000
)
vl_submit(cmd, overwrite = FALSE, logs = "/groups/stark/vloubiere/packages/tests/logs/")

# Example with a BAM file (PRO-Seq data)
cmd <- cmd_demultiplexVBCfile(
  vbcFile = "/groups/stark/projects/PE50_20230401/AAAYLY5HV_1_20230331B_20230401.bam",
  layout = "PAIRED,
  i7 = "GTCCGC",
  i5 = "none,
  i7.column = 14,
  output.prefix = "PE_PROSeq",
  fq.output.folder = "/groups/stark/vloubiere/packages/tests/",
  start.seq = "ATCG",
  trim.length = 10,
  cores = 8,
  head = 40000
)
vl_submit(cmd,
          overwrite = TRUE,
          logs = "/groups/stark/vloubiere/packages/tests/logs/")

}
\seealso{
The underlying Perl script documentation can be accessed using:
\code{system2(system.file("demultiplexing", "vbc_bam_demultiplexing.pl",
             package = "genomicsPipelines"), args = "--help")}
}
