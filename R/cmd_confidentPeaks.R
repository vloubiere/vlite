#' Title
#'
#' @param replicates_peaks_files
#' @param merge_peaks_file
#' @param conf_peaks_output_folder
#' @param Rpath
#'
#' @return
#' @export
#'
#' @examples
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
  if(!dir.exists(conf_peaks_output_folder))
    dir.create(conf_peaks_output_folder, recursive = TRUE, showWarnings = FALSE)

  # Check input files extension ----
  extension <- gsub(".*[.](.*$)", "\\1", merge_peaks_file)

  # Output files paths ----
  peaks.file <- file.path(conf_peaks_output_folder, paste0(output_prefix, ".", extension))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "confident_peaks.R", package = "genomicsPipelines"),
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
