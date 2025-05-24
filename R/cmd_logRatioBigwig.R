#' Call Peaks from bigWig Files
#'
#' @description
#' Generates a shell command to compute the log2 ratio between two bed files using a sliding window.
#' Coverage is computed and normalized for sequencing depth, followed by optional Gaussian smoothing.
#' A pseudocount is added if zeros are present, and log2 ratios are then calculated.
#'
#' @param experiment.bed.file Path to the bed file with experiment reads (must be unique).
#' @param input.bed.file Path to the bed file with input/control reads (must be unique).
#' @param genome BSgenome name.
#' @param gaussian.smoothing If set to TRUE, applies Gaussian smoothing after read-depth normalization. Default= FALSE.
#' @param pseudocount Numeric value to add as a pseudocount if zeros are present after read-depth normalization and
#' optional smoothing. If NULL, the 0.01 percentile of non-zero values is used. Default is NULL.
#' @param bed.subset Optional path to a bed file. If provided, only reads overlapping these regions on the same
#' strand are used. Default= NULL (no filtering).
#' @param bins.width Integer specifying the width using for the sliding window. Default= 100L.
#' @param output.prefix Prefix for output file. If not provided, it is derived from the experiment bed filename.
#' @param output.folder Output directory for the log2 ratio bw file. Default: "db/bw/".
#'
#' @return A data.table with:
#' - `file.type`: Output file label "bw".
#' - `path`: Path to the output file.
#' - `cmd`: Shell command to run.
#' - `job.name`: Default name for the job = "bwLogRatio".
#'
#' @examples
#' cmd <- cmd_logRatioBigwig(
#'   experiment.bed.file = "/data/exp.bed",
#'   input.bed.file = "/data/input.bed",
#' )
#' vl_submit(cmd, execute = FALSE)
#'
#' @export
cmd_logRatioBigwig <- function(experiment.bed.file,
                               input.bed.file,
                               genome,
                               gaussian.smoothing= FALSE,
                               pseudocount= NULL,
                               bed.subset= NULL,
                               bins.width= 100L,
                               output.prefix,
                               output.folder= "db/bw/",
                               Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if files exist to allow wrapping!) ----
  if(length(experiment.bed.file)!=1)
    stop("A unique experiment bed file should be provided.")
  if(length(input.bed.file)!=1)
    stop("A unique input bed file should be provided.")
  if(!is.null(pseudocount) && !is.numeric(pseudocount))
    stop("pseudocount should be numeric.")
  if(any(bins.width %% 1 != 0))
    stop("bins.width should be an integer value.")
  if(missing(output.prefix))
    output.prefix <- gsub(".bw$", "", basename(experiment.bw.file))

  # Output files paths ----
  output.file <- file.path(output.folder, paste0(output.prefix, "_log2_input_ratio.bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "logRatioBigwig.R", package = "vlite"),
    experiment.bed.file, # Experiment bed file
    input.bed.file, # Input bed file
    genome,
    ifelse(is.null(bed.subset), "NULL", bed.subset), # An optional bed file restricting the regions
    gaussian.smoothing,
    ifelse(is.null(pseudocount), "NULL", pseudocount), # Pseudocount to be used
    output.file, # Output file
    bins.width # Bins width
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= output.file,
                    cmd= cmd,
                    job.name= "bwLogRatio")

  # Return ----
  return(cmd)
}
