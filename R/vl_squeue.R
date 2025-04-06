# Squeue to the system
#' @export
vl_squeue <- function()
{
  system("ssh localhost squeue -u vincent.loubiere",
         ignore.stderr = T)
}
