#' Run MAGeCK Analysis for ORFeome Screens
#'
#' @description
#' This function generates shell commands to perform MAGeCK analysis for ORFeome screens.
#' It creates raw and filtered count tables, runs MAGeCK for differential analysis, generates MA plots,
#' and merges results with a master table. The function supports filtering, pseudocount addition, and paired analysis.
#'
#' @param sample.counts A character vector of paths to sample count files.
#' @param input.counts A character vector of paths to input count files.
#' @param screen_name A character string specifying the screen name (e.g., `"ISRE_A549"`). Used to create subfolders in the output directory.
#' @param cdition_name A character string specifying the condition name (e.g., `"IFN_dim"`). Used to create sample output names.
#' @param sample.names A character vector of sample names corresponding to `sample.counts`. Default: extracted from `sample.counts` filenames.
#' @param input.names A character vector of input names corresponding to `input.counts`. Default: extracted from `input.counts` filenames.
#' @param master.table Path to the master table file. Default: `"/groups/stark/pachano/projects/eORFeome_new/Rdata/Master_eORFeome_clean.csv"`.
#' @param output.folder Directory where all output files will be saved. Default: `"db/FC_tables/ORFeome/"`.
#' @param sort A character string specifying the screen type. Options are `"pos"` or `"neg"`. Default: `"pos"`.
#' @param sample.cutoff.FUN A function to filter sample columns. Default: `function(x) sum(x) >= 3`.
#' @param input.cutoff.FUN A function to filter input columns. Default: `function(x) sum(x) >= 0`.
#' @param row.cutoff.FUN A function to filter all columns. Default: `function(x) sum(x) >= 3`.
#' @param pseudocount A numeric value specifying the pseudocount to add to zero values in sample and input columns. Default: `0`.
#' @param paired Logical. Should samples be analyzed as paired? Default: `FALSE`.
#' @param logFC.cutoff A numeric value specifying the logFC cutoff for calling hits in the MA plot. Default: `1`.
#' @param FDR.cutoff A numeric value specifying the FDR cutoff for calling hits in the MA plot. Default: `0.05`.
#' @param Rpath Path to the Rscript executable. Default: `"/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"`.
#'
#' @return A `data.table` with the following columns:
#' - `file.type`: Labels for the output files (e.g., `"raw.counts.table"`, `"filtered.counts.table"`, `"gene.summary"`, `"volcano.plot"`, `"gene.summary.master"`).
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands to run the MAGeCK analysis pipeline.
#'
#' @examples
#' # Example usage
#' cmd <- vl_ORFeome_MAGeCK(
#'   sample.counts = c("sample1_counts.txt", "sample2_counts.txt"),
#'   input.counts = c("input1_counts.txt", "input2_counts.txt"),
#'   screen_name = "ISRE_A549",
#'   cdition_name = "IFN_dim",
#'   sample.names = c("sort_rep1", "sort_rep2"),
#'   input.names = c("input_rep1", "input_rep2"),
#'   output.folder = "db/FC_tables/ORFeome/",
#'   sort = "pos",
#'   pseudocount = 1,
#'   paired = TRUE
#' )
#' print(cmd)
#' vl_submit(cmd, overwrite = FALSE)
#'
#' @export
cmd_MAGeCK_ORFeome <- function(sample.counts,
                               input.counts,
                               screen_name,
                               cdition_name,
                               sample.names= gsub("^(.*?_rep\\d+).*", "\\1", basename(sample.counts)),
                               input.names= gsub("^(.*?_rep\\d+).*", "\\1", basename(input.counts)),
                               master.table= "/groups/stark/pachano/projects/eORFeome_new/Rdata/Master_eORFeome_clean.csv",
                               output.folder= "db/FC_tables/ORFeome/",
                               sort= "pos",
                               sample.cutoff.FUN= function(x) sum(x)>=3,
                               input.cutoff.FUN= function(x) sum(x)>=0,
                               row.cutoff.FUN= function(x) sum(x)>=3,
                               pseudocount= 0,
                               paired= FALSE,
                               logFC.cutoff= 1,
                               FDR.cutoff= 0.05,
                               Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(any(duplicated(sample.names)))
    stop("Some sample.names are duplicated. Make sure that sample.counts files are unique and provide unique sample.names.")
  if(any(duplicated(input.names)))
    stop("Some input.names are duplicated. Make sure that input.counts files are unique and provide unique sample.names.")
  if(!sort %in% c("pos", "neg"))
    stop("sort should either be set to 'pos' or 'neg'.")
  paste0(output_prefix, "_", sort, "FC.txt")

  # Output file names ----
  raw.counts.table <- paste0(output.folder, "/", screen_name, "_", cdition_name, ".", sort, ".raw_counts.txt")
  filtered.counts.table <- gsub(".raw_counts.txt$", ".filtered_counts.txt", raw.counts.table)
  gene.summary <- gsub(".raw_counts.txt$", ".gene_summary.txt", raw.counts.table)
  volcano.plot <- gsub(".raw_counts.txt$", ".volcano_plot.pdf", raw.counts.table)
  gene.summary.master <- gsub("summary.txt$", paste0("summary_master.txt"), gene.summary)

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
    filtered.counts.table # Filtered counts output file
  )
  cmd <- data.table(file.type= c("raw.counts.table", "filtered.counts.table"),
                    path= c(raw.counts.table, filtered.counts.table),
                    cmd= cmd)

  # MAGeCK command ----
  cmd1 <- paste(
    "module load build-env/2020; module load mageck/0.5.9-foss-2018b; mageck test",
    "-k", filtered.counts.table,
    "-t", paste0(unlist(sample.names), collapse = ","),
    "-c", paste0(unlist(input.names), collapse = ","),
    ifelse(paired, "--paired", ""),
    "-n", gsub(".gene_summary.txt", "", gene.summary),
    "--sort-criteria", sort,
    "--remove-zero none"
  )
  cmd <- rbind(cmd,
               data.table(file.type= "gene.summary",
                          path= gene.summary,
                          cmd= cmd1))

  # MA plot and merged table ----
  cmd2 <- paste(
    Rpath,
    system.file("Rscript", "volcano_plots_MAgECK.R", package = "vlite"),
    gene.summary, # Gene summary output file from mageck
    logFC.cutoff, # logFC cutoff for hits
    FDR.cutoff, # FDR cutoff for hits
    volcano.plot # pdf output file
  )
  cmd <- rbind(cmd,
               data.table(file.type= "volcano.plot",
                          path= volcano.plot,
                          cmd= cmd2))

  # Merge master table plot (very fast) ----
  cmd3 <- paste(
    Rpath,
    system.file("Rscript", "merge_gene_summary_to_master_table_ORFeome.R", package = "vlite"),
    gene.summary, # Gene summary output file from mageck \n
    master.table, # 2/ Path to master_table \n
    sort, # sort. can be one of 'pos' or 'neg'. \n
    gene.summary.master # Output file path
  )
  cmd <- rbind(cmd,
               data.table(file.type= "gene.summary.master",
                          path= gene.summary.master,
                          cmd= cmd3))

  # Return ----
  return(cmd)
}
