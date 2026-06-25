#' Filters certain chromosomes from bam file
#'
#' @description
#' Creates shell commands to convert a BAM file to BigWig format. Supports single-end and paired-end sequencing data,
#' with optional fragment extension for paired-end reads.
#'
#' @param bam Path to the input BAM file. Only a single BAM file is allowed.
#' @param seqnames.regexpr A named vector of length 1 containing the regular expression
#' used to filter input bam based on seqnames. The name will be appended to the output.prefix 
#' to create the output basename. Default= c("Dmel"= "^(chr2L|chr2R|chr3L|chr3R|chr4|chrX|chrY)$"); i.e. canonical
#' Dmel chromosomes.
#' @param output.prefix Prefix for the output filtered bam file. If not provided, it is derived from the input BAM filename.
#' @param bam.output.folder Output folder where filtered bam files will be saved.
#' @param cores Number of CPU cores to use. Default= 6.
#'
#' @return A data.table containing the command and output file names.
#'
#' @examples
#'
#' @export
cmd_filterBamChrs <- function(bam,
                              seqnames.regexpr= c("Dmel"= "^(chr2L|chr2R|chr3L|chr3R|chr4|chrX|chrY)$"),
                              output.prefix= NULL,
                              bam.output.folder= "db/bam/",
                              cores= 6)
{
  # Check (!Do not check if bam file exists to allow wrapping!) ----
  if(length(bam)!=1)
    stop("A unique bam file should be provided.")
  if(length(seqnames.regexpr)!=1)
    stop("seqnames.regexpr should be of length one")
  if(is.null(names(seqnames.regexpr)))
    stop("seqnames.regexpr should be named, as the name will be appended to the output prefix to create the output basename")
  if(is.null(output.prefix))
    output.prefix <- gsub(".bam$", "", basename(bam))
  
  # Output files paths ----
  bam.filtered <- file.path(bam.output.folder, paste0(output.prefix, "_", names(seqnames.regexpr), ".filtered.bam"))
  
  # Samtools collapsing command
  filter.cmd <- paste(
    "samtools view -@",
    cores - 1,
    "-h",
    shQuote(bam),
    "| awk",
    shQuote(sprintf('BEGIN{OFS="\t"} /^@/ || $3 ~ /%s/', seqnames.regexpr)),
    "| samtools view -@",
    cores - 1,
    "-b -o",
    shQuote(bam.filtered),
    "-"
  )
  
  # Wrap commands output ----
  cmd <- data.table(file.type= "filtered.bam",
                    path= bam.filtered,
                    cmd= filter.cmd,
                    cores= cores,
                    job.name= "filterBam")
  
  # Return ----
  return(cmd)
}
