#' Squeue to the System
#'
#' This function retrieves the SLURM job queue for a specified user by running the `squeue` command
#' on the SLURM server and reads the output into a data table.
#'
#' @param user A character string specifying the username. Defaults to `"vincent.loubiere"`.
#'
#' @return
#' A data table containing the output of the `squeue` command.
#'
#' @examples
#' \dontrun{
#' # Retrieve the SLURM job queue for the default user
#' vl_squeue()
#'
#' # Retrieve the SLURM job queue for a different user
#' vl_squeue(user = "another.user")
#' }
#'
#' @export
vl_squeue <- function(user = "vincent.loubiere") {
  fread(cmd = paste("ssh localhost squeue -u", user))
}
