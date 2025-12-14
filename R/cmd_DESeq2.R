#' Generate Commands for 'default' DESeq2 analysis
#'
#' @param count.files A comma-separated vector of .txt count files.
#' @param sample.names A comma-separated vector of sample names.
#' @param conditions A comma-separated vector of condition names.
#' @param ctl.conditions comma-separated list of control conditions names.
#' @param norm.counts A comma-separated vector of spike-in/libsize counts that will be used for normalization.
#' Default = NULL (use default method to degine sizeFactors).
#' @param padj.cutoff The p.adjust cutoff to call differentially expressed genes. Default= 0.05.
#' @param log2FC.cutoff The log2FoldChange cutoff to call differentially expressed genes. Default= log2(1.5).
#' @param output.prefix Output prefix to be used for output files.
#' @param dds.output.folder Output folder where DESeq2 dds files ill be same output folder. Default= "db/dds/".
#' @param FC.tables.output.folder Output folder where FC tables will be saved. Default= "db/FC_tables/".
#' @param MAplots.output.folder Output folder where MA plot pdf files will be saved. Default= "pdf/MAplots/".
#' @param Rpath Path to the Rscript binary. Default= "Rscript".
#'
#' @export
cmd_DESeq2 <- function(count.files,
                       sample.names,
                       conditions,
                       ctl.conditions= unique(conditions),
                       norm.counts= NULL,
                       padj.cutoff= 0.05,
                       log2FC.cutoff= log2(1.5),
                       output.prefix,
                       dds.output.folder= "db/dds/",
                       FC.tables.output.folder= "db/FC_tables/",
                       MAplots.output.folder= "pdf/MAplots/",
                       Rpath= "Rscript")
{
  # Check (!Do not check if count.files exist to allow wrapping!) ----
  if(!all(grepl(".txt$", count.files)))
    stop("count.files should all be in .txt format.")
  if(any(!file.exists(count.files)))
    stop("Some file(s) in count.files count not be found.")
  if(uniqueN(lengths(list(count.files, sample.names, conditions)))!=1)
    stop("count.files, sample.names and conditions should all have the same length.")
  if(any(!ctl.conditions %in% conditions))
    stop("All ctl.conditions should exist in conditions")
  if(!is.null(norm.counts)) {
    if(length(norm.counts)!=length(count.files) | !is.numeric(norm.counts))
      stop("norm.counts should be a numeric vector of the same length as count.files.")
    norm.counts <- paste0(norm.counts, collapse = ",")
  }

  # Output files paths ----
  dds.file <- file.path(dds.output.folder, paste0(output.prefix, "_DESeq2.dds"))
  MA.plots <- file.path(MAplots.output.folder, paste0(output.prefix, "_MAplots.pdf"))
  FC.tables <- file.path(FC.tables.output.folder, paste0(output.prefix, "_DESeq2_FC.txt"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "DESeq2_analysis.R", package = "vlite"),
    paste0(count.files, collapse= ","), # A comma-separated list of counts
    paste0(sample.names, collapse= ","), # A comma-separated list of sample names
    paste0(conditions, collapse= ","), # A comma-separated list of condition names
    paste0(ctl.conditions, collapse= ","), # A comma-separated list of control condition names
    padj.cutoff, # p adjust cutoff to call diff genes
    log2FC.cutoff, # log2FC cutoff to call diff genes
    dds.output.folder, # dds output folder
    FC.tables.output.folder, # FC tables output folder
    MAplots.output.folder, # PDF output folder
    output.prefix, # Experiment name
    norm.counts # An optional, comma-separated vector of spike-in/libsize counts for normalization
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("DESeq2.dds", "FC.table", "MA.plot"),
                    path= c(dds.file, FC.tables, MA.plots),
                    cmd= cmd,
                    job.name= "DESeq2")

  # Return ----
  return(cmd)
}
