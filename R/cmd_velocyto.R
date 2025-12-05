#' Generate Cell Ranger Alignment Command
#'
#' @description
#' Constructs a shell command to run velocyto on a Cell Ranger output, producing a .loom file with spliced/unspliced counts from the sample’s BAM and a provided GTF, with resources and modules recorded.
#' Velocyto quantifies spliced and unspliced mRNA from single-cell RNA‑seq to infer RNA velocity—i.e., the short‑term future transcriptional state of cells.
#'
#' @param cellranger.output.folder
#' @param gtf
#' @param cores
#' @param mem
#'
#' @export
cmd_velocyto <- function(cellranger.output.folder,
                         gtf,
                         cores= 8,
                         mem= 64)
{
  # Check ----
  # if(!grepl("outs/$", cellranger.output.folder))
  #   stop("cmd_velocyto: cellranger.output.folder does not end up with '/outs/'")
  if(!grepl(".gtf$", gtf))
    stop("cmd_velocyto: gtf should point to a .gtf file path.")

  # File paths ----
  gz.barcodes <- file.path(cellranger.output.folder, "filtered_feature_bc_matrix", "barcodes.tsv.gz")
  barcodes <- file.path(cellranger.output.folder, "filtered_feature_bc_matrix", "barcodes.tsv")
  bam <- file.path(cellranger.output.folder, "possorted_genome_bam.bam")
  # Output loom file: parent of outs, named after sample
  sample_name <- basename(dirname(cellranger.output.folder))
  loom_dir <- dirname(cellranger.output.folder)
  loom <- file.path(loom_dir, paste0(sample_name, ".loom"))

  # Velocyto command ----
  cmd <- paste(
    if(!file.exists(barcodes)) paste("gunzip -c", gz.barcodes, ">", barcodes, ";") else "",
    "velocyto run",
    "-b", barcodes,
    "-o", loom_dir,
    "-e", sample_name,
    # if (!is.null(repeat_mask)) paste("-m", repeat_mask) else "",
    bam,
    gtf
  )

  # Wrap commands output ----
  cmd <- data.table(file.type= "loom.spliced.unspliced",
                    path= loom,
                    cmd= cmd,
                    cores= cores,
                    mem= mem,
                    job.name= "velocyto",
                    modules= c("build-env/2020",
                               "samtools/1.9-foss-2018b",
                               "velocyto/0.17.17-foss-2018b-python-3.6.6"))

  # Return ----
  return(cmd)
}
