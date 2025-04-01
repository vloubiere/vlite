#' @param fq1
#'
#' @param fq.output.folder
#'
#' @export
cmd_trimProseqAdaptors <- function(fq1,
                                   fq.output.folder= "db/fq/")
{
  # Check ----
  if(length(fq1)!=1)
    stop("A unique fq1 file should be provided.")
  # if(!is.null(fq2) && length(fq2)!=1)
  #   stop("A unique fq2 file should be provided.")
  if(!dir.exists(fq.output.folder))
    dir.create(fq.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output trimmed fq files paths ----
  fq1.trim <- file.path(fq.output.folder, gsub(".fq.gz", "_trimmed.fq.gz", basename(fq1)))

  # Create command ----
  cmd <- paste("cutadapt -a TGGAATTCTCGGGTGCCAAGG",
               "-m", 10, # Remove trimmed reads shorter than 10bp
               "-o", fq1.trim, fq1)

  # Wrap commands output ----
  cmd <- data.table(file.type= "fq1.trim",
                    path= fq1.trim,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
