#' Generate Commands for DESeq2 Analysis of PRO-seq Data
#'
#' @description
#' This function generates shell commands to perform differential expression analysis of PRO-seq data using DESeq2.
#' It supports various normalization methods, including spike-in normalization, and outputs DESeq2 objects, fold-change tables, and MA plots.
#'
#' @param umi.count.tables A character vector of paths to UMI count files in `.txt` format.
#' @param sample.names A character vector of sample names corresponding to the UMI count files.
#' @param conditions A character vector of experimental conditions for each sample.
#' @param ctl.conditions A character vector of control conditions. Default: all unique values in `conditions`.
#' @param ref.genome.stat.files A character vector of paths to reference genome statistics files in `.txt` format.
#' @param spikein.stat.files A character vector of paths to spike-in statistics files in `.txt` format.
#' @param normalization A character string specifying the normalization method. Options are `default`, `libsize`, or `spikeIn`.
#'    Default: `spikeIn`.
#' @param output.prefix A character string specifying the prefix for output files.
#' @param feature A character string specifying the genomic feature being analyzed (e.g., `"exons"`, `"genes"`).
#' @param dds.output.folder Directory for the DESeq2 object output. Default: `db/dds/PROseq/`.
#' @param FC.tables.output.folder Directory for the fold-change tables. Default: `db/FC_tables/PROseq/`.
#' @param MAplots.output.folder Directory for the MA plots. Default: `pdf/MAplots/PROseq/`.
#' @param Rpath Path to the Rscript binary. Default: `/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript`.
#'
#' @return A `data.table` with the following columns:
#' - `file.type`: Labels for the output files (`DESeq2.dds`, `FC.table`, `MA.plot`).
#' - `path`: Paths to the output files.
#' - `cmd`: Shell command to run the DESeq2 analysis pipeline.
#'
#' @examples
#' # Example usage
#' cmd <- vl_PROseq_DESeq2(
#'   umi.count.tables = c("sample1_counts.txt", "sample2_counts.txt"),
#'   sample.names = c("sample1", "sample2"),
#'   conditions = c("treated", "control"),
#'   ctl.conditions = "control",
#'   ref.genome.stat.files = c("sample1_ref_genome_stats.txt", "sample2_ref_genome_stats.txt"),
#'   spikein.stat.files = c("sample1_spikein_genome_stats.txt", "sample2_spikein_genome_stats.txt"),
#'   normalization = "spikeIn",
#'   output.prefix = "experiment1",
#'   feature = "transcripts"
#' )
#' print(cmd)
#'
#' @export
vl_PROseq_DESeq2 <- function(umi.count.tables,
                             sample.names,
                             conditions,
                             ctl.conditions= unique(conditions),
                             ref.genome.stat.files,
                             spikein.stat.files,
                             normalization= "spikeIn",
                             output.prefix,
                             feature,
                             dds.output.folder= "db/dds/PROseq/",
                             FC.tables.output.folder= "db/FC_tables/PROseq/",
                             MAplots.output.folder= "pdf/MAplots/PROseq/",
                             Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(!all(grepl(".txt$", umi.count.files)))
    stop("umi.count.files should all be in .txt format.")
  if(any(!file.exists(umi.count.files)))
    stop("Some file(s) in umi.count.files count not be found.")
  if(uniqueN(lengths(list(umi.count.files, sample.names, conditions, ref.genome.stat.files, spikein.stat.files)))!=1)
    stop("count.files, sample.names, conditions, ref.genome.stat.files and spikein.stat.files should all have the same length.")
  if(!all(grepl(".txt$", ref.genome.stat.files)))
    stop("ref.genome.stat.files should all be in .txt format.")
  if(!all(grepl(".txt$", spikein.stat.files)))
    stop("spikein.stat.files should all be in .txt format.")
  if(any(!ctl.conditions %in% conditions))
    stop("All ctl.conditions should exist in conditions")
  if(length(normalization)!=1)
    stop("normalization should be unique")
  if(!normalization %in% c("default", "libsize", "spikeIn"))
    stop("normalization should be one of 'default', 'libsize', 'spikeIn'")

  # Output files paths ----
  dds.file <- file.path(dds.output.folder, paste0(output.prefix, "_", feature, "_", normalization, "_norm.dds"))
  MA.plots <- file.path(MAplots.output.folder, paste0(output.prefix, "_", feature, "_", normalization, "_MAplots.pdf"))
  FC.tables <- file.path(FC.tables.output.folder, paste0(output.prefix, "_", feature, "_", normalization, "_DESeq2_FC.txt"))

  # DESeq2 command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "DESeq2_PROseq_analysis.R", package = "vlite"),
    paste0(umi.count.tables, collapse = ","), # A comma-separated list of count (ref genome)
    paste0(ref.genome.stat.files, collapse = ","), # A comma-separated list of read statistics (reference genome)
    paste0(spikein.stat.files, collapse = ","), # A comma-separated list of spike-in statistics
    paste0(sample.names, collapse = ","), # A comma-separated list of sample names
    paste0(conditions, collapse = ","), # A comma-separated list of conditions
    paste0(ctl.conditions, collapse = ","), # A comma-separated list of controls
    dds.output.folder, # dds output folder
    FC.tables.output.folder, # FC tables output folder
    MAplots.output.folder, # PDF output folder
    output.prefix, # Experiment
    feature, # feature
    normalization # Normalization method. Possible values are 'default' (DESeq2 default), 'libSize', 'spikeIn', 'combined'.
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("DESeq2.dds", "FC.table", "MA.plot"),
                    path= c(dds.file, FC.tables, MA.plots),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
