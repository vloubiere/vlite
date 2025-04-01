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
cmd_umiToBigwigProseq <- function(umi.counts,
                                  output.prefix= NULL,
                                  bw.output.folder= "db/bw/",
                                  Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript")
{
  # Checks ----
  if(length(umi.counts)!=1)
    stop("A unique umi.counts file should be provided.")
  if(is.null(output.prefix))
    output.prefix <- gsub(".txt$", "", basename(umi.counts))
  if(!dir.exists(bw.output.folder))
    dir.create(bw.output.folder, recursive = TRUE, showWarnings = FALSE)

  # Output file ----
  ps.bw <- file.path(bw.output.folder, paste0(output.prefix, ".ps.bw"))
  ns.bw <- file.path(bw.output.folder, paste0(output.prefix, ".ns.bw"))

  # Command ----
  cmd <- paste(
    Rpath,
    system.file("Rscripts", "umiToBigwigProseq", package = "genomicsPipelines"),
    umi.counts,
    output.prefix
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= c("ps.bw", "ns.bw"),
                    path= c(ps.bw, ns.bw),
                    cmd= cmd)

  # Return ----
  return(cmd)
}
