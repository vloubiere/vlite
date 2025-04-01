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
cmd_bamToBigwig <- function(bam,
                            layout,
                            output.prefix= NULL,
                            extend.PE.fragments= FALSE,
                            extsize= 0,
                            bw.output.folder= "db/bw/",
                            Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(!layout %in% c("SINGLE", "PAIRED"))
    stop("Layout should be one of 'SINGLE' or 'PAIRED'")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  if(extsize<0)
    stop("extsize < 0 not supported")
  if(!dir.exists(bw.output.folder))
    dir.create(bw.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output file ----
  bw.file <- file.path(bw.output.folder, paste0(output.prefix, ".bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "bam_to_bigwig.R", package = "genomicsPipelines"),
    bam,
    layout,
    extend.PE.fragments,
    bw.file,
    extsize
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "bw",
                    path= bw.file,
                    cmd= cmd)

  # Return ----
  return(cmd)
}
