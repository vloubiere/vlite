#' Generate Cell Ranger Aggregate command
#'
#' @description
#'
#' @param molecule_info_paths A character vector of paths to 'molecule_info.h5' files.
#' @param output.prefix Prefix for output files.
#' @param output.folder Directory where output files will be saved. Default= 'db/scRNASeq_10X/aggregated_samples/'.
#' @param cores Number of CPU cores to use. Default= 4.
#' @param mem Memory in GB to allocate. Default= 32.
#'
#' @examples
#'
#' @export
cmd_cellRangerAggr <- function(molecule_info_paths,
                               output.prefix,
                               output.folder= "db/snRNASeq_10X/aggregated_samples/",
                               cores= 4,
                               mem= 32)
{
  # Check
  if(!all(basename(molecule_info_paths)=="molecule_info.h5"))
    stop("The basename of all molecule info paths should be molecule_info.h5")
  molecule_info_paths <- unique(molecule_info_paths)

  # Create csv file ----
  molecule_info_paths <- unique(molecule_info_paths)
  sample_id <- sapply(molecule_info_paths, function(x) rev(unlist(tstrsplit(x, "/")))[3])
  csv <- data.table(sample_id = sample_id,
                    molecule_h5 = molecule_info_paths)

  # Save csv file ----
  csv.file <- file.path(output.folder, paste0(output.prefix, "_samples.csv"))
  fwrite(csv, csv.file, sep= ",")

  # Output file path ----
  output.folder <- file.path(output.folder, output.prefix)
  output.file <- file.path(output.folder, "outs/count/filtered_feature_bc_matrix.h5")

  # Output files paths ----
  cmd <- paste0(
    "cellranger aggr",
    " --id=", output.prefix,
    " --output-dir=", output.folder,
    " --csv=", csv.file
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "h5",
                    path= output.file,
                    cmd= cmd,
                    cores= cores,
                    mem= mem,
                    job.name= "cellRangerAggr",
                    modules= c("build-env/f2022", "cellranger/9.0.0"),
  )

  # Return ----
  return(cmd)
}
