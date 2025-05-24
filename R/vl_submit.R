#' Submit Commands to a Server with Resource Management
#'
#' @description
#' Submits shell commands to a server using resource management (e.g., LSF). Handles job submission,
#' resource allocation, and optional directory creation for output files.
#'
#' @param cmd A data.table containing the commands to execute, with columns:
#'   - `file.type`: Type of output file.
#'   - `path`: Path to the output file.
#'   - `cmd`: Shell command to execute.
#'   - `cores`: (Optional) The number of CPU cores to allocate to the job, specified by the first value.
#'   - `job.name`: (Optional) The name of the job, specified by the first value in this column.
#' @param cores Number of CPU cores to allocate. If not specified, the first value of the optional
#' 'core' column will be used. Otherwise, default= 8.
#' @param mem Memory allocation in GB. Default= 16.
#' @param time Maximum runtime for the job (HH:MM:SS). Default= '08:00:00'.
#' @param job.name Name of the job. If not specified, the first value of the optional 'job.name'
#' column will be used. otherwise, default= "myJob".
#' @param logs Directory for log files. Default: 'db/logs'.
#' @param overwrite If set to TRUE, overwrites existing output files. Default= FALSE.
#' @param execute If set to FALSE, the command is returned and is not submitted to the cluster.
#' Default= TRUE.
#'
#' @examples
#' # Submit a command to the server
#' cmd <- data.table(
#'   file.type = "bam",
#'   path = "output.bam",
#'   cmd = "echo 'Hello world'",
#'   job.name= "test"
#' )
#' vl_submit(cmd)
#'
#' # Return the constructed command without submitting
#' vl_submit(cmd, execute = FALSE)
#'
#' @export
vl_submit <- function(cmd,
                      cores,
                      mem= 16,
                      time= '08:00:00',
                      job.name,
                      logs= "db/logs",
                      overwrite= FALSE,
                      execute= TRUE)
{
  # Checks ----
  if(!is.data.table(cmd) || !all(c("file.type", "path", "cmd") %in% names(cmd)))
    stop("cmd should be a data.table containing columns 'file.type', 'path' and 'cmd'.")

  # Default job name and cores ----
  if(missing(job.name)) {
    job.name <- if("job.name" %in% names(cmd))
      cmd$job.name[1] else
        "myJob"
  }
  if(missing(cores)) {
    cores <- if("cores" %in% names(cmd))
      cmd$cores[1] else
        8
  }

  # Check existing files ----
  if(!overwrite) {
    cmd <- cmd[!file.exists(path)]
  }

  # Check if any command should be executed ----
  if(nrow(cmd))
  {
    # Create missing directories ----
    dirs <- unique(dirname(cmd$path))
    if(any(!dir.exists(dirs)) && execute) {
      sapply(dirs[!dir.exists(dirs)], dir.create, recursive = TRUE, showWarnings = FALSE)
    }

    # Create final command ----
    cmd <- paste(c("module load build-env/2020",
                   "module load cutadapt/1.18-foss-2018b-python-2.7.15",
                   "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                   "module load samtools/1.9-foss-2018b",
                   "module load bowtie/1.2.2-foss-2018b",
                   "module load bowtie2/2.3.4.2-foss-2018b",
                   "module load macs2/2.1.2.1-foss-2018b-python-2.7.15",
                   "module load mageck/0.5.9-foss-2018b",
                   "module load macs2/2.1.2.1-foss-2018b-python-2.7.15",
                   unique(cmd$cmd)), collapse = "; ")
    if(execute) {
      # Submit ----
      jobID <- bsub(cmd= cmd,
                    cores= cores,
                    mem = mem,
                    time = time,
                    name = job.name,
                    logs= logs)
      print(jobID)
    } else {
      return(cmd)
    }
  } else
    message("All output files exist and overwrite= FALSE. No command submitted.")
}
