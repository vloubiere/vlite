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
#' @param pdf.output.folder Directory for the pdf containing diagnostic plots. Default= "db/peaks/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#'
#' @return A data.table with:
#' - `file.type`: Output file label ("confident.peaks").
#' - `path`: Path to the confident peaks file.
#' - `cmd`: Shell command to run the confident peaks pipeline.
#' - `job.name`: Default name for the job = "confPeaks".
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
                               pdf.output.folder= "db/peaks/",
                               Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Check (!Do not check if replicates.peaks.files exist to allow wrapping!) ----
  if(!all(grepl(".narrowPeak$|broadPeak$", c(replicates.peaks.files, merge.peaks.file))))
    stop("replicates.peaks.files and merge.peaks.file should be in .narrowPeak or .broadPeaks format.")
  if(!all(grepl(".narrowPeak$|broadPeak$", merge.peaks.file)))
    stop("merge.peaks.file should be in .narrowPeak or .broadPeaks format.")

  # Check input files extension ----
  extension <- gsub(".*[.](.*$)", "\\1", merge.peaks.file)

  # Output files paths ----
  peaks.file <- file.path(conf.peaks.output.folder, paste0(output.prefix, ".", extension))
  pdf.file <- file.path(pdf.output.folder, paste0(output.prefix, ".pdf"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscript", "confident_peaks.R", package = "vlite"),
    paste0(replicates.peaks.files, collapse= ","), # Replicates
    merge.peaks.file, # Merged reads
    peaks.file, # Diagnostic pdf
    pdf.file # Diagnostic pdf
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("confident.peaks", "confident.peaks.pdf"),
                    path= c(peaks.file, pdf.file),
                    cmd= cmd,
                    job.name= "confPeaks")

  # Return ----
  cmd$mem <- 4
  cmd$cores <- 1
  return(cmd)
}
