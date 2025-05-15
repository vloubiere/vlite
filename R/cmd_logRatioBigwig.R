#' Call Peaks from bigWig Files
#'
#' @description
#' Generates a shell command to compute log2 ratio between two bed files, using a sliding window.
#'
#' @param experiment.bed.file Path to a bed file containing experiment reads. Must be unique.
#' @param input.bed.file Path to a bed file containing input/control reads. Must be unique.
#' @param genome A BS genome name.
#' @param pseudocount Pseudocount to use and avoid 0s. If set to NULL, it will be set to the 0.01 percentile
#' of non-0 values. Default= NULL.
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
    ifelse(is.null(pseudocount), "NULL", pseudocount), # Pseudocount to be used
    output.file, # Output file
    bins.width # Bins width
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= output.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
