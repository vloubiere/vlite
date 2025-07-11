% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmd_umiToBigwigProseq.R
\name{cmd_umiToBigwigProseq}
\alias{cmd_umiToBigwigProseq}
\title{Convert UMI Counts to BigWig Format for PRO-Seq Data}
\usage{
cmd_umiToBigwigProseq(
  umi.counts,
  output.prefix = NULL,
  bw.output.folder = "db/bw/",
  Rpath = "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"
)
}
\arguments{
\item{umi.counts}{Path to the input UMI counts file. Only a single file is allowed.}

\item{output.prefix}{Prefix for the output BigWig files. If not provided, it is derived from the
input UMI counts filename.}

\item{bw.output.folder}{Directory for the BigWig files. Default= "db/bw/".}

\item{Rpath}{Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".}
}
\value{
A data.table with:
- `file.type`: Output file labels ("ps.bw", "ns.bw").
- `path`: Paths to the positive strand and negative strand BigWig files.
- `cmd`: Shell command to run the UMI to BigWig conversion.
- `job.name`: Default name for the job = "umiToBwProseq".
}
\description{
Creates shell commands to convert UMI counts from PRO-Seq data into positive strand and negative
strand, CPM=-normalized bigwig files.
}
\examples{
# Convert UMI counts to BigWig files for PRO-Seq data
cmd <- cmd_umiToBigwigProseq(
  umi.counts = "/data/umi_counts/sample_UMI_counts.txt"
)
vl_submit(cmd, execute= FALSE)

}
