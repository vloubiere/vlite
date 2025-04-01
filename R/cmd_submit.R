#' Title
#'
#' @param cmd
#' @param cores
#' @param mem
#' @param time
#' @param job.name
#' @param logs
#' @param overwrite
#' @param Rpath
#' @param execute Should the command be submitted to the server? Default= TRUE
#'
#' @return
#' @export
#'
#' @examples
cmd_submit <- function(cmd,
                       cores= 8,
                       mem= 16,
                       time= '08:00:00',
                       job.name= "myJob",
                       logs= "db/logs",
                       overwrite= FALSE,
                       Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                       execute= TRUE)
{
  # Checks ----
  if(!is.data.table(cmd) || !all(c("file.type", "path", "cmd") %in% names(cmd)))
    stop("cmd should be a data.table containing columns 'file.type', 'path' and 'cmd'.")

  # Check which commands should be executed ----
  if(!overwrite) {
    cmd <- cmd[!file.exists(path)]
  }

  if(nrow(cmd))
  {
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
      bsub(cmd= cmd,
           cores= cores,
           mem = mem,
           time = time,
           name = job.name,
           logs= logs)
    } else {
      return(cmd)
    }

  } else
    message("All output files already existed! No command submitted ;). Consider overwrite= T if convenient.")
}
