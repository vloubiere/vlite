#' Process RNA-Seq Data
#'
#' @description
#' Implements a pipeline for processing RNA-Seq data, including:
#' 1. Adapter trimming
#' 2. Genome alignment
#' 3. Gene counting
#' 4. BigWig track generation
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths.
#' @param fq2 For paired-end data, a character vector of .fq (or .fq.gz) file paths matching fq1 files. Default= NULL.
#' @param output.prefix Prefix for output files.
#' @param genome Reference genome identifier (e.g., "mm10", "hg38").
#' @param genome.idx Path to the Rsubread genome index. If NULL, derived from genome. Default= NULL.
#' @param gtf Path to the GTF annotation file. If genome is provided, corresponding protein-coding / mRNA genes will be used. Default= NULL.
#' @param GTF.attrType.extra Additional GTF attribute to include in the output (e.g., gene symbol).
#' @param fq.output.folder Directory for trimmed FASTQ files. Default= "db/fq/RNASeq/".
#' @param bam.output.folder Directory for aligned BAM files. Default= "db/bam/RNASeq/".
#' @param alignment.stats.output.folder Directory for alignment statistics. Default= "db/alignment_stats/RNASeq/".
#' @param counts.stats.output.folder Directory for gene count statistics. Default= "db/stats/RNASeq/".
#' @param counts.output.folder Directory for gene counts. Default= "db/counts/RNASeq/".
#' @param bw.output.folder Directory for BigWig files. Default= "db/bw/RNASeq/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of CPU cores to use. Default= 8.
#'
#' @return A data.table with:
#' - `file.types`: Types of output files.
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands for each step in the pipeline.
#'
#' @examples
#' # Process paired-end RNA-Seq data
#' cmd <- rnaseqProcessing(
#'   fq1 = c("sample1_R1.fq.gz", "sample1_R1_rep2.fq.gz"),
#'   fq2 = c("sample1_R2.fq.gz", "sample1_R2_rep2.fq.gz"),
#'   output.prefix = "sample1",
#'   genome = "hg38",
#'   gtf = "/data/annotations/hg38.gtf",
#'   GTF.attrType.extra = "gene_name",
#'   cores = 8
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Then, differential analysis should be performed using ?cmd_DESeq2()
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{cmd_trimIlluminaAdaptors}} for adapter trimming
#'   \item \code{\link{cmd_alignRsubread}} for alignment
#'   \item \code{\link{cmd_countRsubread}} for gene counting
#'   \item \code{\link{cmd_bamToBigwig}} for BigWig track generation
#' }
#'
#' @export
rnaseqProcessing <- function(fq1,
                             fq2= NULL,
                             output.prefix,
                             genome,
                             genome.idx= NULL,
                             gtf= NULL,
                             GTF.attrType.extra,
                             fq.output.folder= "db/fq/RNASeq/",
                             bam.output.folder= "db/bam/RNASeq/",
                             alignment.stats.output.folder= "db/alignment_stats/RNASeq/",
                             counts.stats.output.folder= "db/stats/RNASeq/",
                             counts.output.folder= "db/counts/RNASeq/",
                             bw.output.folder= "db/bw/RNASeq/",
                             Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                             cores= 8)
{
  # Trimming illumina adaptors ----
  cmd <- cmd_trimIlluminaAdaptors(fq1= fq1,
                                  fq2= fq2,
                                  fq.output.folder= fq.output.folder)

  # * If several fq1/fq2 files provided, they will be merged during alignment ----
  fq1.trim <- cmd[file.type=="fq1.trim", path]
  if(!is.null(fq2))
    fq2.trim <- cmd[file.type=="fq2.trim", path]

  # Alignment ----
  align.cmd <- cmd_alignRnaRsubread(fq1= fq1.trim,
                                    fq2= fq2.trim,
                                    output.prefix= output.prefix,
                                    genome= genome,
                                    genome.idx= genome.idx,
                                    bam.output.folder= bam.output.folder,
                                    alignment.stats.output.folder = alignment.stats.output.folder,
                                    Rpath= Rpath)
  cmd <- rbind(cmd, align.cmd)

  # Gene counts ----
  count.cmd <- cmd_countRsubread(bam= cmd[file.type=="bam", path],
                                 layout= ifelse(is.null(fq2), "SINGLE", "PAIRED"),
                                 output.prefix = NULL, # bam file basename
                                 genome= genome,
                                 gtf= gtf,
                                 GTF.attrType.extra= GTF.attrType.extra,
                                 counts.stats.output.folder= counts.stats.output.folder,
                                 counts.output.folder= counts.output.folder,
                                 Rpath= Rpath)
  cmd <- rbind(cmd, count.cmd)

  # bw tracks ----
  bw.cmd <- cmd_bamToBigwig(bam = align.cmd[file.type=="bam", path],
                            layout = ifelse(is.null(fq2), "SINGLE", "PAIRED"),
                            output.prefix = NULL, # bam file basename
                            extend.PE.fragments = FALSE,
                            extsize = 0,
                            bw.output.folder = bw.output.folder,
                            Rpath = Rpath)
  cmd <- rbind(cmd, bw.cmd)

  # DESeq2 should be performed separately, as controls cannot be fetched here ----

  # Return ----
  return(cmd)
}
