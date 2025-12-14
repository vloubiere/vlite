#' Run MAGeCK Analysis for ORFeome Screens
#'
#' @description
#' This function generates shell commands to perform MAGeCK analysis for ORFeome screens.
#' It creates raw and filtered count tables, runs MAGeCK for differential analysis, generates MA plots,
#' and merges results with a master table. The function supports filtering, pseudocount addition, and paired analysis.
#'
#' @param sample.counts A character vector of paths to sample count files.
#' @param sample.names For each sample.counts file, the name of the corresponding sample. By default, it will be inferred from basename(sample.counts).
#' @param input.counts A character vector of paths to input count files.
#' @param input.names For each input.counts file, the name of the corresponding sample. By default, it will be inferred from basename(input.counts).
#' @param sample.cutoff.FUN A function to filter sample rows. Default: `function(x) sum(x) >= 3`.
#' @param input.cutoff.FUN A function to filter input rows. Default: `function(x) sum(x) >= 0`.
#' @param row.cutoff.FUN A function to filter rows using both sample and input columns. Default: `function(x) sum(x) >= 3`.
#' @param pseudocount A numeric value specifying the pseudocount to add to zero values in sample and input columns. Default= 0.
#' @param sort A character string specifying whether positive ('pos') or negative ('neg') enrichments are expected. Default= 'pos'.
#' @param paired Logical. Should samples be analyzed as paired? Default= FALSE.
#' @param logFC.cutoff A numeric value specifying the logFC cutoff for calling hits. Default= log2(1.5).
#' @param FDR.cutoff A numeric value specifying the FDR cutoff for calling hits. Default= 0.05.
#' @param master.table.rds Path to a .rds file containing the master table that will be merged to MAGeCK output.
#' Should at least contain the ORF ids in a column named 'id'.
#' @param output.prefix Output files prefix (e.g. 'screen1_sample1_vs_input').
#' @param FC.output.folder Directory where final FC tables should be saved. Default= 'db/FC_tables/ORFeome/'.
#' @param MAGeCK.output.folder Directory where MAGeCK intermediate files should be saved. Default= 'db/FC_tables/ORFeome/'.
#' @param Rpath Path to the Rscript executable. Default: "Rscript".
#'
#' @return A data.table containing columns:
#' - `file.type`: Labels for the output files (i.e., "raw.counts.table", "filtered.counts", "gene.summary", "volcano.plot", "gene.summary.master").
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands to run the MAGeCK analysis pipeline.
#' - `job.name`: Default name for the job = "MAGECK".
#'
#' @examples
#'
#' @export
cmd_MAGeCK_ORFeome <- function(sample.counts,
                               sample.names= gsub("^(.*?_rep\\d+).*", "\\1", basename(sample.counts)),
                               input.counts,
                               input.names= gsub("^(.*?_rep\\d+).*", "\\1", basename(input.counts)),
                               sample.cutoff.FUN= function(x) sum(x)>=3,
                               input.cutoff.FUN= function(x) sum(x)>=0,
                               row.cutoff.FUN= function(x) sum(x)>=3,
                               pseudocount= 0,
                               sort= "pos",
                               paired= FALSE,
                               logFC.cutoff= log2(1.5),
                               FDR.cutoff= 0.05,
                               master.table.rds,
                               output.prefix,
                               FC.output.folder= "db/FC_tables/ORFeome/",
                               MAGeCK.output.folder= "db/FC_tables/ORFeome/MAGeCK/",
                               Rpath= "Rscript")
{
  # Checks (!Do not check if input files exist to allow wrapping) ----
  if(any(duplicated(sample.counts)))
    stop("Some sample.counts paths are duplicated.")
  if(any(duplicated(sample.names)))
    stop("Some sample.names are duplicated.")
  if(any(duplicated(input.counts)))
    stop("Some input.counts paths are duplicated.")
  if(any(duplicated(input.names)))
    stop("Some input.names are duplicated.")
  if(!sort %in% c("pos", "neg"))
    stop("sort should either be set to 'pos' or 'neg'.")

  # Output file names ----
  raw.counts.table <- file.path(MAGeCK.output.folder, paste0(output.prefix, ".", sort, ".raw_counts.txt"))
  filtered.counts <- file.path(MAGeCK.output.folder, paste0(output.prefix, ".", sort, ".filtered_counts.txt"))
  normalized.counts <- file.path(MAGeCK.output.folder, paste0(output.prefix, ".", sort, ".normalized.txt"))
  gene.summary <- file.path(MAGeCK.output.folder, paste0(output.prefix, ".", sort, ".gene_summary.txt"))
  FC.table <- file.path(FC.output.folder, paste0(output.prefix, ".", sort, "_FC_MAGeCK.txt"))
  volcano.plot <- file.path(FC.output.folder, paste0(output.prefix, ".", sort, ".volcano_plot.pdf"))

  # Deparse functions ----
  sample_fun_string <- paste(deparse(sample.cutoff.FUN), collapse = " ")
  input_fun_string <- paste(deparse(input.cutoff.FUN), collapse = " ")
  row_fun_string <- paste(deparse(row.cutoff.FUN), collapse = " ")

  # Command to generate raw and filtered counts table ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "compute_MAGECK_count_tables_ORFeome.R", package = "vlite"),
    paste0(sample.counts, collapse = ","), # A comma-separated list of sample count files
    paste0(input.counts, collapse = ","), # A comma-separated list of input count files
    paste0(sample.names, collapse = ","), # A comma-separated list of sample names (mathcing the list of count files)
    paste0(input.names, collapse = ","), # A comma-separated list of input names (mathcing the list of count files)
    shQuote(sample_fun_string), # Function to be applied to input columns for filtering
    shQuote(input_fun_string), # Function to be applied to sample columns for filtering
    shQuote(row_fun_string), # Function to be applied to all columns for filtering
    pseudocount, # Pseudocount to be added to input and sample columns (which will be further normalized for sequencing depth)
    raw.counts.table, # Raw counts output file
    filtered.counts # Filtered counts output file
  )
  cmd <- data.table(file.type= c("raw.counts.table", "filtered.counts"),
                    path= c(raw.counts.table, filtered.counts),
                    cmd= cmd)

  # MAGeCK command ----
  cmd1 <- paste(
    "module load build-env/2020; module load mageck/0.5.9-foss-2018b; mageck test",
    "-k", filtered.counts,
    "-t", paste0(unlist(sample.names), collapse = ","),
    "-c", paste0(unlist(input.names), collapse = ","),
    "--sort-criteria", sort,
    "--norm-method median",
    "--normcounts-to-file",
    ifelse(paired, "--paired", ""),
    "-n", gsub(".gene_summary.txt", "", gene.summary),
    "--remove-zero none"
  )
  cmd <- rbind(cmd,
               data.table(file.type= c("raw.counts", "normalized.counts", "gene.summary"),
                          path= c(raw.counts.table, normalized.counts, gene.summary),
                          cmd= cmd1))

  # Merge gene summary to Master table ----
  cmd2 <- paste(
    Rpath,
    system.file("Rscript", "merge_gene_summary_to_master_table_ORFeome.R", package = "vlite"),
    gene.summary, # 1/ Gene summary output file from MAGeCK (.txt).
    sort, # 2/ sort. can be one of 'pos' or 'neg'.
    logFC.cutoff, # 3/ logFC cutoff to define hits.
    FDR.cutoff, # 4/ FDR cutoff to define hits.
    master.table.rds, # 5/ Path to .rds master_table
    FC.table # 6/ Output file path (.txt).
  )
  cmd <- rbind(cmd,
               data.table(file.type= "FC.table",
                          path= FC.table,
                          cmd= cmd2))

  # MA plot and merged table ----
  cmd3 <- paste(
    Rpath,
    system.file("Rscript", "volcano_plots_MAgECK.R", package = "vlite"),
    FC.table, # Final FC table
    volcano.plot # pdf output file
  )
  cmd <- rbind(cmd,
               data.table(file.type= "volcano.plot",
                          path= volcano.plot,
                          cmd= cmd3))

  # Return ----
  cmd$job.name <- "MAGECK"
  return(cmd)
}
