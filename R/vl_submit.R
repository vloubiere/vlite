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
#' @param cores Number of CPU cores to allocate. If the cmd object contains a 'cores' column, it will override this argument. Default= 8.
#' @param mem Memory allocation in GB. If the cmd object contains a 'mem' column, it will override this argument. Default= 32.
#' @param time Maximum runtime for the job (HH:MM:SS). If the cmd object contains a 'time' column, it will override this argument. Default= '08:00:00'.
#' @param job.name Name of the job. If the cmd object contains a 'job.name' column, it will override this argument. Default= "myJob".
#' @param logs Directory for log files. If the cmd object contains a 'logs' column, it will override this argument. Default = 'db/logs'.
#' @param modules A character vector of modules to be loaded before executing the command. If the cmd object contains a 'modules' column,
#' it will override this argument. Default= c("build-env/f2022", "mamba/24.3.0-0"). By default, the conda modules used by the pipelines are loaded.
#' @param conda If specified, the conda environment will be activated after loading the modules. If the cmd object contains a 'conda' column,
#' it will override this argument. Default= "/groups/stark/conda_envs/.conda/envs/vlite_pipelines/".
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
                      cores= 8,
                      mem= 32,
                      time= '08:00:00',
                      job.name= "myJob",
                      logs= "db/logs/",
                      modules= c("build-env/f2022", "mamba/24.3.0-0"),
                      conda= "/groups/stark/conda_envs/.conda/envs/vlite_pipelines/",
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
  if('cores' %in% names(cmd))
    cores <- cmd$cores[1]
  if('mem' %in% names(cmd))
    mem <- cmd$mem[1]
  if('time' %in% names(cmd))
    time <- cmd$time[1]
  if('job.name' %in% names(cmd))
    job.name <- cmd$job.name[1]
  if('logs' %in% names(cmd))
    logs <- cmd$logs[1]
  if('modules' %in% names(cmd))
    modules <- unlist(cmd$modules)
  if('conda' %in% names(cmd))
    conda <- cmd$conda[1]

  # ensure no empty or duplicated modules ----
  modules <- unique(trimws(unlist(modules)))
  modules <- modules[!is.na(modules) & nzchar(modules)]

  # Check if output files exist ----
  if(!overwrite)
    cmd <- cmd[is.na(path) | !file.exists(path)]

  # If any command was provided ----
  if(nrow(cmd))
  {
    # Create missing directories ----
    valid_paths <- !is.na(cmd$path) & nzchar(cmd$path)
    dirs <- unique(dirname(cmd$path[valid_paths]))
    if(any(!dir.exists(dirs)) && create.output.dirs && execute) {
      sapply(dirs[!dir.exists(dirs)], dir.create, recursive = TRUE, showWarnings = FALSE)
    }

    # Create final command ----
    final_cmd <- paste(unique(cmd$cmd), collapse = "; ")
    if(!is.null(conda)) {
      if(!any(grepl("mamba", modules)))
        warning("A conda env was specified but no mamba module was loaded!")
      final_cmd <- paste0("source activate ", conda, "; ", final_cmd)
    }
    if(length(modules))
      final_cmd <- paste(c(paste0("ml ", modules), final_cmd), collapse = "; ")

    # Submit ----
    if(execute) {
      jobID <- bsub(
        cmd= final_cmd,
        cores= cores,
        mem = mem,
        time = time,
        name = job.name,
        logs= logs
      )
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
