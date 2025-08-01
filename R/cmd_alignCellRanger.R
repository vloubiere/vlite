#' Generate Cell Ranger Alignment Command
#'
#' @description
#' Constructs a shell command to run 10x Genomics Cell Ranger for single-cell RNA-seq alignment and quantification, with resource specification and output management.
#'
#' @param fq.prefix Prefix shared between input fastq files: '/path/to/fastq/prefix'.
#' cellranger will only search for file names matching 'prefix_S*_L.*'.
#' @param index Path to the Cell Ranger reference transcriptome. e.g.: '/groups/stark/vloubiere/genomes/Drosophila_melanogaster/cellranger_dm6_TG/index/'
#' @param GEX.multiome Do the fastq files correspond to the GEX part of a multiome experiment? Default= F
#' @param output.prefix Prefix for output files.
#' @param output.folder Directory where output files will be saved. Default= 'db/scRNASeq_10X/'.
#' @param cores Number of CPU cores to use. Default= 10.
#' @param mem Memory in GB to allocate. Default= 64.
#'
#' @examples
#' cmd_alignCellRanger(
#'   fq.prefix = "/data/fastq/sample1",
#'   index = "/data/genomes/cellranger_index",
#'   output.prefix = "sample1",
#'   output.folder = "db/scRNASeq_10X/"
#' )
#'
#' @export
cmd_alignCellRanger <- function(fq.prefix,
                                index,
                                GEX.multiome= FALSE,
                                output.prefix,
                                output.folder= "db/scRNASeq_10X/",
                                cores= 10,
                                mem= 64)
{
  # Check (!Do not check if fq1 or fq2 files exist to allow wrapping!) ----
  fq.prefix <- unique(fq.prefix)
  if(length(fq.prefix)>1)
    stop("Error: unique(fq.prefix)>1.
    All fq files should be in the same folder and share the same prefix before 'S*_L.*'.
         If necessary, rename them consistently.")
  # Check number of fastq files
  fq.dir <- dirname(fq.prefix)
  fq.basename <- basename(fq.prefix)
  check <- list.files(fq.dir, fq.basename)
  if(length(check)<2 | !all(grepl(".fq$|.fastq$|.fq.gz$|.fastq.gz$", check)))
    stop("Less than two fq files were found.") else
      print(paste0(length(check), " fastq files found!"))

  # Output files paths ----
  bam <- file.path(output.folder, paste0(output.prefix, ".bam"))

  # Output files paths ----
  cmd <- paste0(
    "module load build-env/f2022;",
    "module load cellranger/9.0.0;",
    "cd ", output.folder, ";",
    " cellranger count --id=", output.prefix,
    " --transcriptome=", index,
    " --fastqs=", normalizePath(fq.dir),
    " --sample=", fq.basename,
    " --create-bam=true",
    " --localcores=", cores,
    " --localmem=", mem
  )

  # To analyze only the GEX part of a multiome dataset ----
  if(GEX.multiome)
    cmd <- paste0(cmd, " --chemistry=ARC-v1")

  # Wrap commands output ----
  cmd <- data.table(file.type= "bam",
                    path= bam,
                    cmd= cmd,
                    cores= cores,
                    mem= mem,
                    job.name= "cellRanger")

  # Return ----
  return(cmd)
}
