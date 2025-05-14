#' Call Peaks from bigWig Files
#'
#' @description
#' Generates a shell command to generate log2 ratio bw file to compared an experiment be track
#'  to an input/control track.
#'
#' @param experiment.bw.file Path to the experiment bigWig (.bw) file. Must be unique.
#' @param input.bw.file Path to input/control bigWig (.bw) file. Must be unique.
#' @param bed Optional path to a bed file to which the ratio calulation will be restricted.
#' If missing, all chromosomes present in experiment and input will be considered.
#' @param pseudocount Pseudocount to use and avoid 0s. If left empty, it will be set to 1% percentile of non-0 values.
#' @param output.prefix Prefix for output file.
#' @param output.folder Output directory for the log2 ratio bw file. Default: "db/bw/".
#'
#' @return A data.table with:
#' - `file.type`: Output file label "bw".
#' - `path`: Path to the output file.
#' - `cmd`: Shell command to run.
#'
#' @examples
#' cmd <- cmd_logRatioBigwig(
#'   experiment.bw.file = "/data/exp.bw",
#'   input.bw.file = "/data/input.bw",
#' )
#' vl_submit(cmd, execute = FALSE)
#'
#' @export
cmd_logRatioBigwig <- function(experiment.bw.file,
                               input.bw.file,
                               bed,
                               pseudocount,
                               output.prefix,
                               output.folder= "db/bw/",
                               Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if files exist to allow wrapping!) ----
  if(length(experiment.bw.file)!=1)
    stop("A unique experiment bw file should be provided.")
  if(length(input.bw.file)!=1)
    stop("A unique input bw file should be provided.")
  if(!missing(pseudocount) && !is.numeric(pseudocount))
    stop("pseudocount should be numeric.")

  # Output files paths ----
  output.file <- file.path(output.folder, paste0(output.prefix, ".bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "logRatioBigwig.R", package = "vlite"),
    experiment.bw.file, # Experiment bw file
    input.bw.file, # Input bw file
    ifelse(missing(bed), "NULL", bed), # An option bed file containing the regions for which log ratio will be computed
    ifelse(missing(pseudocount), "NULL", pseudocount), # The pseudocount to be used
    output.file # Output file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= output.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
