#' Title
#'
#' @param bam
#' @param layout
#' @param output.prefix
#' @param extend.PE.fragments
#' @param bw.output.folder
#' @param extsize
#' @param Rpath
#'
#' @return
#' @export
#'
#' @examples
cmd_bedToBigwig <- function(bed,
                            genome,
                            output.prefix= NULL,
                            bw.output.folder= "db/bw/",
                            Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(bed)!=1)
    stop("A unique bed file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bed$", "", basename(bed))
  if(!dir.exists(bw.output.folder))
    dir.create(bw.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output file ----
  bw.file <- file.path(bw.output.folder, paste0(output.prefix, "_", genome, ".bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "bed_to_bigwig.R", package = "genomicsPipelines"),
    bed,
    genome,
    bw.file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= bw.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
