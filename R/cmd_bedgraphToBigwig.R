#' Title
#'
#' @param bdg
#' @param output.prefix
#' @param genome
#' @param Rpath
#'
#' @return
#' @export
#'
#' @examples
cmd_bedgraphToBigwig <- function(bdg,
                                 output.prefix,
                                 bw.output.folder= "db/bw/",
                                 genome,
                                 Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(bdg)!=1)
    stop("Only one bedgraph file should be provided.")
  if(!dir.exists(bw.output.folder))
    dir.create(bw.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output file ----
  bw.file <- file.path(bw.output.folder, paste0(output.prefix, "_", genome, ".bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "bedgraph_to_bigwig.R", package = "genomicsPipelines"),
    bdg,
    bw.file,
    genome
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= bw.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
