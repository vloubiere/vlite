# Scancel all jobs except Rstudio
#' @export
vl_scancel <- function()
{
  .c <- fread(cmd= "ssh localhost squeue -u vincent.loubiere")
  .c <- .c[!(NAME %in% c("[RStudio", "jupyter_"))]
  if(nrow(.c))
    system(paste(c("ssh localhost scancel", .c$JOBID), collapse= " "),
           ignore.stderr = T)
}
