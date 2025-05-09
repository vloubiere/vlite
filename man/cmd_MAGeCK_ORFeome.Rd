% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmd_MAGECK_ORFeome.R
\name{cmd_MAGeCK_ORFeome}
\alias{cmd_MAGeCK_ORFeome}
\title{Run MAGeCK Analysis for ORFeome Screens}
\usage{
cmd_MAGeCK_ORFeome(
  sample.counts,
  input.counts,
  screen_name,
  cdition_name,
  sample.names = gsub("^(.*?_rep\\\\d+).*", "\\\\1", basename(sample.counts)),
  input.names = gsub("^(.*?_rep\\\\d+).*", "\\\\1", basename(input.counts)),
  master.table =
    "/groups/stark/pachano/projects/eORFeome_new/Rdata/Master_eORFeome_clean.csv",
  output.folder = "db/FC_tables/ORFeome/",
  sort = "pos",
  sample.cutoff.FUN = function(x) sum(x) >= 3,
  input.cutoff.FUN = function(x) sum(x) >= 0,
  row.cutoff.FUN = function(x) sum(x) >= 3,
  pseudocount = 0,
  paired = FALSE,
  logFC.cutoff = 1,
  FDR.cutoff = 0.05,
  Rpath = "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"
)
}
\arguments{
\item{sample.counts}{A character vector of paths to sample count files.}

\item{input.counts}{A character vector of paths to input count files.}

\item{screen_name}{A character string specifying the screen name (e.g., `"ISRE_A549"`). Used to create subfolders in the output directory.}

\item{cdition_name}{A character string specifying the condition name (e.g., `"IFN_dim"`). Used to create sample output names.}

\item{sample.names}{A character vector of sample names corresponding to `sample.counts`. Default: extracted from `sample.counts` filenames.}

\item{input.names}{A character vector of input names corresponding to `input.counts`. Default: extracted from `input.counts` filenames.}

\item{master.table}{Path to the master table file. Default: `"/groups/stark/pachano/projects/eORFeome_new/Rdata/Master_eORFeome_clean.csv"`.}

\item{output.folder}{Directory where all output files will be saved. Default: `"db/FC_tables/ORFeome/"`.}

\item{sort}{A character string specifying the screen type. Options are `"pos"` or `"neg"`. Default: `"pos"`.}

\item{sample.cutoff.FUN}{A function to filter sample columns. Default: `function(x) sum(x) >= 3`.}

\item{input.cutoff.FUN}{A function to filter input columns. Default: `function(x) sum(x) >= 0`.}

\item{row.cutoff.FUN}{A function to filter all columns. Default: `function(x) sum(x) >= 3`.}

\item{pseudocount}{A numeric value specifying the pseudocount to add to zero values in sample and input columns. Default: `0`.}

\item{paired}{Logical. Should samples be analyzed as paired? Default: `FALSE`.}

\item{logFC.cutoff}{A numeric value specifying the logFC cutoff for calling hits in the MA plot. Default: `1`.}

\item{FDR.cutoff}{A numeric value specifying the FDR cutoff for calling hits in the MA plot. Default: `0.05`.}

\item{Rpath}{Path to the Rscript executable. Default: `"/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"`.}
}
\value{
A `data.table` with the following columns:
- `file.type`: Labels for the output files (e.g., `"raw.counts.table"`, `"filtered.counts.table"`, `"gene.summary"`, `"volcano.plot"`, `"gene.summary.master"`).
- `path`: Paths to the output files.
- `cmd`: Shell commands to run the MAGeCK analysis pipeline.
}
\description{
This function generates shell commands to perform MAGeCK analysis for ORFeome screens.
It creates raw and filtered count tables, runs MAGeCK for differential analysis, generates MA plots,
and merges results with a master table. The function supports filtering, pseudocount addition, and paired analysis.
}
\examples{
# Example usage
cmd <- vl_ORFeome_MAGeCK(
  sample.counts = c("sample1_counts.txt", "sample2_counts.txt"),
  input.counts = c("input1_counts.txt", "input2_counts.txt"),
  screen_name = "ISRE_A549",
  cdition_name = "IFN_dim",
  sample.names = c("sort_rep1", "sort_rep2"),
  input.names = c("input_rep1", "input_rep2"),
  output.folder = "db/FC_tables/ORFeome/",
  sort = "pos",
  pseudocount = 1,
  paired = TRUE
)
print(cmd)
vl_submit(cmd, overwrite = FALSE)

}
