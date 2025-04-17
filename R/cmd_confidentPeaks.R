#' Generate Commands for Identifying Confident Peaks
#'
#' @description
#' Creates shell commands to identify confident peaks across replicates using merged peak files.
#' Outputs a file containing confident peaks.
#'
#' @param replicates.peaks.files Vector of paths to replicate peak files in .narrowPeak or .broadPeak format.
#' @param merge.peaks.file Path to the merged peak file in .narrowPeak or .broadPeak format.
#' @param output.prefix Prefix for the output file.
#' @param conf.peaks.output.folder Directory for the confident peaks file. Default= "db/peaks/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file label ("confident.peaks").
#' - `path`: Path to the confident peaks file.
#' - `cmd`: Shell command to run the confident peaks pipeline.
#'
#' @examples
#' # Identify confident peaks across replicates
#' cmd <- cmd_confidentPeaks(
#'   replicates.peaks.files = c("/data/peaks/rep1.narrowPeak", "/data/peaks/rep2.narrowPeak"),
#'   merge.peaks.file = "/data/peaks/merged.narrowPeak",
#'   output.prefix = "sample1"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_confidentPeaks <- function(replicates.peaks.files,
                               merge.peaks.file,
                               output.prefix,
                               conf.peaks.output.folder= "db/peaks/",
                               Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if replicates.peaks.files exist to allow wrapping!) ----
  if(!all(grepl(".narrowPeak$|broadPeak$", replicates.peaks.files)))
    stop("replicates.peaks.files should be in .narrowPeak or .broadPeaks format.")
  if(!all(grepl(".narrowPeak$|broadPeak$", merge.peaks.file)))
    stop("merge.peaks.file should be in .narrowPeak or .broadPeaks format.")

  # Check input files extension ----
  extension <- gsub(".*[.](.*$)", "\\1", merge.peaks.file)

  # Output files paths ----
  peaks.file <- file.path(conf.peaks.output.folder, paste0(output.prefix, ".", extension))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "confident_peaks.R", package = "vlite"),
    replicates.peaks.files, # Replicates
    merge.peaks.file, # Merged reads
    peaks.file # Output file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "confident.peaks",
                    path= peaks.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
