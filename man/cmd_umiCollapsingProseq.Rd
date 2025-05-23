% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmd_umiCollapsingProseq.R
\name{cmd_umiCollapsingProseq}
\alias{cmd_umiCollapsingProseq}
\title{Generate Commands to collapse UMIs from a bam file containing PRO-Seq reads}
\usage{
cmd_umiCollapsingProseq(
  bam,
  output.prefix = NULL,
  counts.output.folder = "db/counts/",
  stats.output.folder = "db/counts_statistics/",
  Rpath = "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"
)
}
\arguments{
\item{bam}{Path to the input BAM file. Only a single BAM file is allowed.}

\item{output.prefix}{Prefix for the output files. If not provided, it is derived from the input BAM filename.}

\item{counts.output.folder}{Directory for the UMI counts file. Default= "db/counts/".}

\item{stats.output.folder}{Directory for the UMI statistics file. Default= "db/counts_statistics/".}

\item{Rpath}{Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".}
}
\value{
A data.table with:
- `file.type`: Output file labels ("umi.counts", "umi.stats").
- `path`: Paths to the output files.
- `cmd`: Shell command to run the UMI counting pipeline.
- `job.name`: Default name for the job = "umiCollProseq".
}
\description{
Creates shell commands to count unique molecular identifiers (UMIs) from a PRO-Seq BAM file.
Outputs a UMI counts file and a statistics file.
}
\examples{
# Count UMIs in a PRO-Seq BAM file
cmd <- cmd_umiCollapsingProseq(
  bam = "/data/bam/sample.bam"
)
vl_submit(cmd, execute= FALSE)

}
