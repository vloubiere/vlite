#' Title
#'
#' @param seurat_object
#' @param loom.file
#' @param cdition
#' @param output.folder
#'
#' @return
#' @export
#'
#' @examples
cmd_filterLoom <- function(seurat_object,
                           input.loom.file= NULL,
                           cdition,
                           output.folder= "db/velocity/")
{
  # Check ----
  stopifnot(grepl(".rds$", seurat_object))
  stopifnot(grepl(".loom$", input.loom.file))

  # Output file ----
  output.file <- file.path(output.folder, cdition, paste0(cdition, "_filtered.loom"))

  # Make command ----
  cmd <- paste(
    "Rscript /groups/stark/vloubiere/projects/scRNASeq_ED/git_scRNASeq/functions/filter_loom_file.R",
    seurat_object,
    input.loom.file,
    cdition,
    output.folder
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "loom",
                    path= output.file,
                    cmd= cmd,
                    cores= 8,
                    mem= 32,
                    job.name= "filterLoom",
                    modules= c("build-env/f2022", "r-bundle-bioconductor/3.19-foss-2023b-r-4.4.1"))

  # Return ----
  return(cmd)
}
