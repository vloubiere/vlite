#' Compute candidate regulons per gene
#'
#' @description
#'
#' @param genes.file An .rds file containing SCT-normalized counts for the genes to be tested (dependent variables).
#' @param TFs.file An .rds file containing SCT-normalized counts for the predictor TFs.
#' @param cell.labels.file An .rds file containing a data.table containing "cluster" and "cdition" columns.
#' @param N.iterations The number of iterations to be performed. Default= 100L.
#' @param cell.frac The fractions of cells to be used per iteration (random sampling). Default= .3.
#' @param stability.cutoff The fraction of time a TF should receive a coefficient to be considered as a candidate
#' activator/repressor. Default= 0.6.
#' @param method Should be one of 'LASSO' (quick) or 'ELASTIC' (less affected by colinearity between predictors but slower).
#' @param output.prefix The prefix of the output file. If missing, the basename of the genes_file will be used.
#' @param cores Number of cores to be used. Default= 20.
#' @param mem Memory dedicated to the job. Default= 40.
#' @param time Time limit for the job. Default= "8:00:00".
#'
#' @export
cmd_infer_candidate_regulons <- function(
    genes.file,
    TFs.file,
    cell.labels.file,
    N.iterations= 100L,
    cell.frac= .3,
    stability.cutoff= .6,
    method= "ELASTIC",
    output.prefix,
    output.folder= "db/candidate_regulons/",
    cores= 20,
    mem= 40,
    time= "8:00:00"
)
{
  # Check ----
  if(missing(output.prefix))
    output.prefix <- gsub(".rds$", "", basename(genes.file))
  stopifnot(grepl(".rds$", genes.file))
  stopifnot(grepl(".rds$", TFs.file))
  stopifnot(grepl(".rds$", cell.labels.file))
  stopifnot(method %in% c("LASSO", "ELASTIC"))

  # File paths ----
  output.file <- file.path(
    output.folder,
    paste0(output.prefix, "_candidate_regulons_", method, ".rds")
  )

  # Primary commands ----
  cmd <- paste(
    "Rscript /groups/stark/vloubiere/vlite/inst/Rscript/infer_candidate_regulons.R",
    genes.file,
    TFs.file,
    cell.labels.file,
    method,
    N.iterations, # N iterations
    cell.frac, # Fraction of random sampled cells
    stability.cutoff, # Fraction of times a candidate TF must receive a coefficient
    output.prefix,
    output.folder
  )

  # Command object ----
  cmd <- data.table(
    cmd= cmd,
    file.type= "cand_regulon",
    path= output.file,
    cores= cores,
    mem= mem,
    time= time,
    modules= list(c("build-env/f2022", "r-bundle-bioconductor/3.19-foss-2023b-r-4.4.1"))
  )

  # Return ----
  return(cmd)
}
