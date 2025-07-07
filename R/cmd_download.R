#' Download files
#'
#' @param url URL to download
#' @param output.name The output file name (default= keep url basename).
#' @param submit If set to TRUE, the command will directly be submitted. If set to FALSE,
#' a 'cmd' object compatible with ?vl_submit will be returned. Default= TRUE.
#' @param output.folder Path where the output file will be saved.
#' @param output.file.type The type of output file (fq...). Default= NA.
#' @param gzip Should the output file be gzipped? Default= FALSE
#'
#' @return A `data.table` with:
#' - `file.type`: Output file labels.
#' - `path`: Paths to trimmed FASTQ files.
#' - `cmd`: Shell command to run Trim Galore.
#' - `job.name`: Default name for the job = "download".
#'
#' @examples
#' cmd <- cmd_download(url= "url/to/data")
#' vl_submit(cmd, execute= FALSE)
#'
#' @export
cmd_download <- function(url,
                         output.name= basename(url),
                         submit= TRUE,
                         output.folder,
                         output.file.type= NA,
                         gzip= FALSE)
{
  # Checks ----
  if(length(url)!=1)
    stop("A unique URL should be provided.")
  if(!dir.exists(dirname(output.folder)))
    dir.create(dirname(output.folder), recursive = TRUE, showWarnings = FALSE)

  # Download command ----
  output.file.path <- file.path(output.folder, output.name)
  cmd <- paste("wget", url, "-O", output.file.path)

  # Automatically gzip the output files ----
  if(gzip) {
    cmd <- paste(cmd, "; gzip", output.file.path)
    output.file.path <- paste0(output.file.path, ".gz")
  }

  # Wrap commands output ----
  cmd <- data.table(file.type= output.file.type,
                    path= output.file.path,
                    cmd= cmd,
                    job.name= "download")

  # Submit or return command ----
  if(submit) {
    vl_submit(cmd, cores = 1, mem = 1, job.name = "dl", logs = "~/")
  } else {
    # Return
    return(cmd)
  }
}
