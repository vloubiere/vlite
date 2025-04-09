#' Cancel SLURM Jobs Except Specified Ones
#'
#' Cancels all SLURM jobs for a user, except those with names matching the specified patterns.
#'
#' @param user Username for which jobs should be canceled. Defaults to `vincent.loubiere`.
#' @param ignore.name Names of the jobs that should not be cancelled. Defaults to c('[RStudio', 'jupyter_').
#'
#' @return
#' No return value. Executes the `scancel` command to cancel jobs.
#'
#' @examples
#' \dontrun{
#' vl_scancel()  # Cancel jobs for the default user
#' vl_scancel(user = "another.user", ignore.name = c("my_job"))
#' }
#'
#' @export
vl_scancel <- function(user= "vincent.loubiere", ignore.name= c("[RStudio", "jupyter_"))
{
  # Squeue
  .c <- vl_squeue(user= user)
  # Remove names to ignore
  .c <- .c[!(NAME %in% ignore.name)]
  # Cancel jobs
  if(nrow(.c))
    system(paste(c("ssh localhost scancel", .c$JOBID), collapse= " "),
           ignore.stderr = T)
}
