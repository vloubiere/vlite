#' Generate Commands for 'default' DESeq2 analysis
#'
#' @export
cmd_DESeq2 <- function(count.files,
                       sample.names,
                       conditions,
                       ctl.conditions= unique(conditions),
                       output.prefix,
                       dds.output.folder= "db/dds/RNASeq/",
                       FC.tables.output.folder= "db/FC_tables/RNASeq/",
                       MAplots.output.folder= "pdf/MAplots/RNASeq/",
                       Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(!all(grepl(".txt$", count.files)))
    stop("count.files should all be in .txt format.")
  if(any(!file.exists(count.files)))
    stop("Some file(s) in count.files count not be found.")
  if(uniqueN(lengths(list(count.files, sample.names, conditions)))!=1)
    stop("count.files, sample.names and conditions should all have the same length.")
  if(any(!ctl.conditions %in% conditions))
    stop("All ctl.conditions should exist in conditions")

  # Output files paths ----
  dds.file <- file.path(dds.output.folder, paste0(output.prefix, ".dds"))
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
    dds.output.folder, # dds output folder
    FC.tables.output.folder, # FC tables output folder
    MAplots.output.folder, # PDF output folder
    output.prefix # Experiment name
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("DESeq2.dds", "FC.table", "MA.plot"),
                    path= c(dds.file, FC.tables, MA.plots),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
