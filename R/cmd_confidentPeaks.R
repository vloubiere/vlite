#' Generate Commands for Identifying Confident Peaks
#'
#' @description
#' Creates shell commands to identify confident peaks across replicates using merged peak files.
#' Outputs a file containing confident peaks.
#'
#' @param replicates_peaks_files Vector of paths to replicate peak files in `.narrowPeak` or `.broadPeak` format.
#' @param merge_peaks_file Path to the merged peak file in `.narrowPeak` or `.broadPeak` format.
#' @param output_prefix Prefix for the output file.
#' @param conf_peaks_output_folder Directory for the confident peaks file. Default: `"db/peaks/"`.
#' @param Rpath Path to the Rscript binary. Default: `"/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"`.
#'
#' @return A `data.table` with:
#' - `file.type`: Output file label (`"confident.peaks"`).
#' - `path`: Path to the confident peaks file.
#' - `cmd`: Shell command to run the confident peaks pipeline.
#'
#' @examples
#' # Identify confident peaks across replicates
#' cmd <- cmd_confidentPeaks(
#'   replicates_peaks_files = c("/data/peaks/rep1.narrowPeak", "/data/peaks/rep2.narrowPeak"),
#'   merge_peaks_file = "/data/peaks/merged.narrowPeak",
#'   output_prefix = "sample1"
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_confidentPeaks <- function(replicates_peaks_files,
                               merge_peaks_file,
                               output_prefix,
                               conf_peaks_output_folder= "db/peaks/",
                               Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(!all(grepl(".narrowPeak$|broadPeak$", replicates_peaks_files)))
    stop("Replicates_peaks_files should be in .narrowPeak or .broadPeaks format.")
  if(!all(grepl(".narrowPeak$|broadPeak$", merge_peaks_file)))
    stop("merge_peaks_file should be in .narrowPeak or .broadPeaks format.")

  # Check input files extension ----
  extension <- gsub(".*[.](.*$)", "\\1", merge_peaks_file)

  # Output files paths ----
  peaks.file <- file.path(conf_peaks_output_folder, paste0(output_prefix, ".", extension))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "confident_peaks.R", package = "vlite"),
    replicates_peaks_files, # Replicates
    merge_peaks_file, # Merged reads
    peaks.file # Output file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "confident.peaks",
                    path= peaks.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
