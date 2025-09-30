#' Extracts the first 10,000 of a VBC tar file
#'
#' @description
#' Generates the call to a perl script to extract the head of VBC tar files for debugging purposes
#'
#' @param vbcTarFile Path to a vbc tar.gz file.
#' @param fq.output.folder Output folder where fq files should be saved. Default= "/scratch-cbe/users/vincent.loubiere/<tar_file_basename>/".
#' @param nreads Number of reads that should be processed (for testing purposes). Default= 10,000L.
#' @param cores Number of CPU cores to use for samtools processing (when using BAM input). Default= 8L.
#' @param mem Memory to use for the job (Go). Default= 16.
#' @param overwrite If output files already exist, should they be overwritten? Default= FALSE.
#'
#' @return A data.table with three columns:
#' \itemize{
#'   \item `file.type`: Labels for output files (e.g.: "fq1", "fq2").
#'   \item `path`: Full paths to the output files.
#'   \item `cmd`: Shell command(s) to generate the output files.
#'   \item `cores`: The number of CPU cores to use.
#'   \item `job.name`: Default name for the job = "demultVBC".
#' }
#'
#' @details
#' @export
cmd_checkVBCfile <- function(vbcTarFile,
                             fq.output.folder= paste0("/scratch-cbe/users/vincent.loubiere/", gsub(".tar.gz", "", basename(vbcTarFile))),
                             nreads= 10000L,
                             cores= 8,
                             mem= 16,
                             overwrite= FALSE)
{
  # Checks ----
  options(scipen= 999)
  # Input file
  if(!grepl(".tar.gz$", vbcTarFile))
    stop("vbcTarFile should be in .bam or .tar.gz format.")
  if(length(vbcTarFile)!=1 || !file.exists(vbcTarFile))
    stop("vbcTarFile does not exist.")
  # Head
  if(nreads %% 1!=0 | nreads<=0)
    stop("nreads should be a round number > 0.")

  # Check if output folder already contains files ----
  paths <- if(dir.exists(fq.output.folder))
    list.files(fq.output.folder, ".fastq$", full.names = T) else
      character(0)

  # Main function ----
  if(length(paths)==0 || overwrite) {

    # Generate command ----
    perl.script <- system.file("perl", "vbc_tar_extract_head.pl", package = "vlite")
    cmd <- paste("cd", fq.output.folder,
                 "perl", perl.script,
                 paste0("--tar ", vbcTarFile),
                 paste0("--reads ", nreads),
                 paste0("--threads ", cores))

    # Create a dummy file (so that vl_submit creates the output folder)
    if(length(paths)==0)
      paths <- file.path(fq.output.folder, "output.files.fastq")

    # Command object ----
    cmd <- data.table(file.type= "fq",
                      path= paths,
                      cmd= cmd,
                      cores= cores,
                      job.name= "checkVBC")

    # Submit ----
    vl_submit(cmd,
              cores = cores,
              mem= mem,
              logs = file.path(fq.output.folder, "logs"),
              overwrite = TRUE)
    print(paste("Output files will be saved in", fq.output.folder))
  } else {
    print(paste("Overwrite was set to FALSE and fastq files already exist in", shQuote(fq.output.folder)))
    print("---> No command was submitted")
  }
}
