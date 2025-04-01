#' Submit Jobs to VBC Server Using bsub
#'
#' @description
#' A wrapper function to submit jobs to the VBC server using bsub/GridEngine.
#' The function handles job submission parameters and can either execute the job
#' or return the formatted command string.
#'
#' @param cmd character, The command to be executed on the server
#' @param cores integer, Number of CPU cores to request (default: 6)
#' @param mem numeric, Memory requirement in GB (default: 12)
#' @param time character, Maximum wall time in format 'HH:MM:SS' (default: '08:00:00')
#' @param name character, Job name prefix (default: "vl")
#' @param logs character, Path to directory for storing log files (default: "db/logs/")
#' @param wdir character, Working directory for job execution (default: current working directory)
#' @param execute logical, If TRUE submits the job, if FALSE returns the formatted command (default: TRUE)
#'
#' @return If execute=TRUE, returns a character string with the job ID.
#'         If execute=FALSE, returns the formatted bsub command string.
#'
#' @details
#' The function creates a bsub command using the provided parameters and either:
#' 1. Submits the job through SSH to localhost (when execute=TRUE)
#' 2. Returns the formatted command string (when execute=FALSE)
#'
#' The function automatically unsets SBATCH environment variables before submission
#' and creates a temporary script file that is removed after job submission.
#'
#' @examples
#' \dontrun{
#' # Submit a simple job
#' bsub("echo 'Hello World'", cores=1, mem=1)
#'
#' # Get the formatted command without executing
#' bsub("echo 'Hello World'", execute=FALSE)
#'
#' # Submit a job with custom parameters
#' bsub(cmd="Rscript analysis.R",
#'      cores=8,
#'      mem=16,
#'      time='12:00:00',
#'      name="analysis",
#'      logs="path/to/logs")
#' }
#'
#' @export
bsub <- function(cmd,
                 cores= 6L,
                 mem= 12,
                 time= '08:00:00',
                 name= "vl",
                 logs= "db/logs/",
                 wdir= getwd(),
                 execute= T)
{
  # Checks ----
  if(!dir.exists(logs))
    dir.create(logs, recursive = TRUE, showWarnings = FALSE)

  # Create tmp folder to store sh files ----
  tmp.folder <- file.path(normalizePath(logs), "tmp")
  if(!dir.exists(tmp.folder))
    dir.create(tmp.folder, recursive = TRUE, showWarnings = FALSE)

  # Assemble command ----
  cmd <- paste0("cd ", wdir,
                "; /groups/stark/software-all/shell/bsub_gridengine",
                " -C ", cores,
                " -n ", name,
                " -m ", mem,
                " -o ", logs,
                " -e ", logs,
                " -T ", time,
                " \"", cmd, "\"")

  # Submit ----
  if(execute) {
    Sys.unsetenv("SBATCH_RESERVATION")
    Sys.unsetenv("SBATCH_WCKEY")
    # Write command in temp file
    sh.file <- tempfile(tmpdir = tmp.folder,
                        fileext = ".sh")
    writeLines(cmd, sh.file)
    # Execute using ssh
    job_ID <- system(paste0("ssh localhost sh ", sh.file),
                     intern = T,
                     ignore.stderr = TRUE)
    # Remove temp file
    file.remove(sh.file)
    return(gsub(".* (.*$)", "JobID: \\1", job_ID[2]))
  } else {
    return(cmd)
  }
}
