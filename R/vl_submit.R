#' Submit Commands to a Server with Resource Management
#'
#' @description
#' Used in vlite pipeline to submit shell commands to the server Handles job submission, resource allocation, and optional directory
#' creation for output files.
#'
#' @param cmd A data.table containing the commands to execute, which will all be submitted as a single job. Intended columns are:
#'   - `file.type`: Type of output file.
#'   - `path`: Path to the output file.
#'   - `cmd`: Shell command to execute.
#'   - `cores`: (Optional) The number of CPU cores to allocate to the job (only the value from the first row will be used).
#'   - `mem`: (Optional) Memory allocation in GB (only the value from the first row will be used).
#'   - `time`: (Optional) Maximum runtime for the job (HH:MM:SS, only the value from the first row will be used).
#'   - `job.name`: (Optional) The name of the job, only the value from the first row will be used.
#'   - `logs`: (Optional) Directory for log files, only the value from the first row will be used
#'   - `modules`: (Optional) A vector or list of modules to be loaded before executing the command.
#' @param cores Number of CPU cores to allocate. If not specified, the first value of cmd$cores column will be used. Default= 8.
#' @param mem Memory allocation in GB. If not specified, the first value of cmd$mem column will be used. Default= 32.
#' @param time Maximum runtime for the job (HH:MM:SS). If not specified, the first value of cmd$time column will be used. Default= '08:00:00'.
#' @param job.name Name of the job. If not specified, the first value of cmd$job.name column will be used. Default= "myJob".
#' @param logs Directory for log files. If not specified, the first value of cmd$logs column will be used. Default = 'db/logs'.
#' @param modules A character vector of modules to be loaded before executing the command. If not specified, the modules listed in cmd$modules will be used.
#' By default, the most used modules for genomics are loaded.
#' @param overwrite If set to TRUE, overwrites existing output files. Default= FALSE.
#' @param execute If set to FALSE, the command is returned and is not submitted to the cluster. Default= TRUE.
#' @param create.output.dirs Should missing output directories be created before executing the command? Default= TRUE.
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
                      modules,
                      overwrite= FALSE,
                      execute= TRUE,
                      create.output.dirs= TRUE)
{
  # Checks ----
  if(!is.data.table(cmd) || !all(c("file.type", "path", "cmd") %in% names(cmd)))
    stop("cmd should be a data.table containing columns 'file.type', 'path' and 'cmd'.")
  if(any(is.na(cmd$path)))
    warning("Some rows have path = NA; they will be skipped for dir creation and existence checks.")

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
      as.character(cmd$logs)[1] else
        "db/logs"
  }
  if(missing(modules)) {
    modules <- if("modules" %in% names(cmd))
      unlist(cmd$modules) else c(
        "build-env/2020",
        "pigz/2.4-gcccore-7.3.0 ",
        "cutadapt/1.18-foss-2018b-python-2.7.15",
        "trim_galore/0.6.0-foss-2018b-python-2.7.15",
        "samtools/1.9-foss-2018b",
        "bowtie/1.2.2-foss-2018b",
        "bowtie2/2.3.4.2-foss-2018b",
        "macs2/2.1.2.1-foss-2018b-python-2.7.15",
        "mageck/0.5.9-foss-2018b"
      )
  }

  # ensure no empty or duplicated modules ----
  modules <- unique(trimws(unlist(modules)))
  modules <- modules[!is.na(modules) & nzchar(modules)]

  # Check if output files exist ----
  if(!overwrite)
    cmd <- cmd[is.na(path) | !file.exists(path)]

  # Check if any command should be executed ----
  if(nrow(cmd))
  {
    # Create missing directories ----
    valid_paths <- !is.na(cmd$path) & nzchar(cmd$path)
    dirs <- unique(dirname(cmd$path[valid_paths]))
    if(any(!dir.exists(dirs)) && create.output.dirs && execute) {
      sapply(dirs[!dir.exists(dirs)], dir.create, recursive = TRUE, showWarnings = FALSE)
    }

    # Create final command ----
    final_cmd <- if(length(modules))
      paste(c(paste0("ml ", modules), unique(cmd$cmd)), collapse = "; ") else
        paste(unique(cmd$cmd), collapse = "; ")

    # Submit ----
    if(execute) {
      jobID <- bsub(cmd= final_cmd,
                    cores= cores,
                    mem = mem,
                    time = time,
                    name = job.name,
                    logs= logs)
      print(jobID)
      invisible(jobID)
    } else {
      return(final_cmd)
    }
  } else {
    message("All output files exist and overwrite= FALSE. No command submitted.")
    invisible(character(0))
  }
}
