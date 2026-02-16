#' Title
#'
#' @param cmd 
#' @param logs 
#' @param name 
#' @param wdir 
#'
#' @returns
#' @export
#'
#' @examples
submitCav <- function(
    cmd,
    conda= "vl_conda",
    logs= "db/logs/",
    name= "vl",
    create.output.dirs= TRUE,
    overwrite= FALSE,
    execute= TRUE,
    wdir= getwd()
)
{
  # Check if output files exist ----
  if(!overwrite)
    cmd <- cmd[is.na(path) | !file.exists(path)]
  
  # If output files are missing ----
  if(nrow(cmd)) {
    
    # Create missing directories ----
    # Output files
    valid_paths <- !is.na(cmd$path) & nzchar(cmd$path)
    dirs <- unique(dirname(cmd$path[valid_paths]))
    if(any(!dir.exists(dirs)) && create.output.dirs && execute) {
      sapply(dirs[!dir.exists(dirs)], dir.create, recursive = TRUE, showWarnings = FALSE)
    }
    # logs
    base::dir.create(logs, showWarnings = FALSE, recursive = TRUE)
    
    # Script and log file output names ----
    tmp <- base::tempfile(pattern = paste0(name, "_"), tmpdir = logs)
    out <- paste0(tmp, ".out")
    err <- paste0(tmp, ".err")
    res <- list(script = tmp, out = out, err = err)
    
    # Assemble final command ----
    # shebang
    cmds <- c(
      "#!/usr/bin/env bash",
      "set -eo pipefail",
      paste0("cd ", base::shQuote(wdir))
    )
    # Load conda env
    if (!is.null(conda)) {
      cmds <- c(
        cmds,
        "source /home/michael.szalay/anaconda3/etc/profile.d/conda.sh",
        "set +u",
        paste("conda activate", conda),
        "set -u"
      )
    }
    
    # Save script file ----
    base::writeLines(paste(c(cmds, unique(cmd$cmd)), collapse = "\n\n"), tmp)
    
    # Execute ----
    if(execute) {
      status <- base::system2(
        command = "bash",
        args = c(tmp),
        stdout = out,
        stderr = err
      )
      # Save status ----
      res <- c(res, list(status))
    }
    
    # Return status ----
    invisible(res) 
    
  } else
    print("All output files already exist.")
}