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
#'   - `mem`: (Optional) Memory allocation in GB.
#'   - `time`: (Optional) Maximum runtime for the job (HH:MM:SS).
#'   - `job.name`: (Optional) The name of the job, specified by the first value in this column.
#'   - `logs`: (Optional) Directory for log files.
#' @param cores Number of CPU cores to allocate. If not specified, the first value of cmd$cores column will be used. Default= 8.
#' @param mem Memory allocation in GB. If not specified, the first value of cmd$mem column will be used. Default= 32.
#' @param time Maximum runtime for the job (HH:MM:SS). If not specified, the first value of cmd$time column will be used. Default= '08:00:00'.
#' @param job.name Name of the job. If not specified, the first value of cmd$job.name column will be used. Default= "myJob".
#' @param logs Directory for log files. If not specified, the first value of cmd$job.name column will be used. Default= 'db/logs'.
#' @param modules A character vector of modules to be loaded before executing the command. By default, the most used modules are loaded.
#' @param overwrite If set to TRUE, overwrites existing output files. Default= FALSE.
#' @param execute If set to FALSE, the command is returned and is not submitted to the cluster. Default= TRUE.
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
                      mem,
                      time,
                      job.name,
                      logs,
                      modules= c(
                        "module load build-env/2020",
                        "module load cutadapt/1.18-foss-2018b-python-2.7.15",
                        "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                        "module load samtools/1.9-foss-2018b",
                        "module load bowtie/1.2.2-foss-2018b",
                        "module load bowtie2/2.3.4.2-foss-2018b",
                        "module load macs2/2.1.2.1-foss-2018b-python-2.7.15",
                        "module load mageck/0.5.9-foss-2018b",
                        "module load macs2/2.1.2.1-foss-2018b-python-2.7.15"
                      ),
                      overwrite= FALSE,
                      execute= TRUE)
{
  # Checks ----
  if(!is.data.table(cmd) || !all(c("file.type", "path", "cmd") %in% names(cmd)))
    stop("cmd should be a data.table containing columns 'file.type', 'path' and 'cmd'.")

  # Default values ----
  if(missing(cores)) {
    cores <- if("cores" %in% names(cmd))
      cmd$cores[1] else
        8
  }
  if(missing(mem)) {
    mem <- if("mem" %in% names(cmd))
      cmd$mem[1] else
        32
  }
  if(missing(time)) {
    time <- if("time" %in% names(cmd))
      cmd$time[1] else
        '08:00:00'
  }
  if(missing(job.name)) {
    job.name <- if("job.name" %in% names(cmd))
      cmd$job.name[1] else
        "myJob"
  }
  if(missing(logs)) {
    logs <- if("logs" %in% names(cmd))
      cmd$logs[1] else
        "db/logs"
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
    cmd <- paste(c(modules, unique(cmd$cmd)), collapse = "; ")
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
