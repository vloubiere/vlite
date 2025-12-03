#' Process PRO-Seq Sequencing Data
#'
#' @description
#' Implements a pipeline for processing PRO-Seq sequencing data, including:
#' 1. Adapter trimming
#' 2. Removal of reads aligning to ncRNAs (decoy).
#' 3. Reference genome alignment
#' 4. Spike-in genome alignment
#' 5. Genome alignment (reference and spike-in)
#' 6. UMI counting
#' 7. BigWig track generation
#'
#' @param fq1 A character vector of .fq (or .fq.gz) file paths. Note that the UMI sequences should be appended to the readIDs
#' (see ?cmd_demultiplexVBCfile()), and that read 2 are not used.
#' @param output.prefix Prefix for output files.
#' @param ref.genome Reference genome identifier (e.g., "mm10", "hg38").
#' @param ref.genome.idx Path to the Bowtie1 index for the reference genome. If NULL, derived from ref.genome. Default= NULL.
#' @param spike.genome Spike-in genome identifier.
#' @param spike.genome.idx Path to the Bowtie1 index for the spike-in genome. Default= NULL.
#' @param ref.gtf Path to the GTF annotation file for the reference genome. Default= NULL.
#' @param fq.output.folder Directory for trimmed FASTQ files. Default= "db/fq/PROSeq/".
#' @param bam.output.folder Directory for aligned BAM files. Default= "db/bam/PROSeq/".
#' @param alignment.stats.output.folder Directory for alignment statistics. Default= "db/alignment_stats/PROSeq/".
#' @param counts.output.folder Directory for UMI counts. Default= "db/counts/PROSeq/".
#' @param counts.stats.output.folder Directory for UMI count statistics. Default= "db/stats/PROSeq/".
#' @param bw.output.folder Directory for BigWig files. Default= "db/bw/PROSeq/".
#' @param Rpath Path to the Rscript binary. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of CPU cores to use. Default= 8.
#' @param mem Memory to use. Default= 64.
#' @param time Max time for job. Default= '2-00:00:00'.
#'
#' @return A data.table with:
#' - `file.types`: Types of output files.
#' - `path`: Paths to the output files.
#' - `cmd`: Shell commands for each step in the pipeline.
#' - `cores`: The number of CPU cores to use.
#' - `job.name`: Default name for the job = paste0("PRO_", output.prefix).
#'
#' @examples
#' # Process PRO-Seq sequencing data
#' cmd <- proseqProcessing(
#'   fq1 = c("sample1_R1.fq.gz", "sample1_R1_rep2.fq.gz"),
#'   output.prefix = "sample1",
#'   ref.genome = "mm10",
#'   spike.genome = "dm6",
#'   ref.gtf = "/data/annotations/mm10.gtf",
#'   cores = 8
#' )
#' vl_submit(cmd, execute= FALSE)
#'
#' # Next, reads should be counted for genomic features of interest using ?cmd_countPROseqReads()
#' # Then, diff analysis should be performed using ?cmd_DESeq2_PROseq()
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{cmd_trimProseqAdaptors}} for adapter trimming
#'   \item \code{\link{cmd_alignBowtie}} for alignment
#'   \item \code{\link{cmd_exractUnalignedReads}} for extracting unaligned reads
#'   \item \code{\link{cmd_umiToBigwigProseq}} for UMI counting and BigWig generation
#' }
#'
#' @export
proseqProcessing_ncRNAdecoy <- function(
    fq1,
    output.prefix,
    ref.genome,
    ref.genome.ncRNA.idx= NULL,
    ref.genome.idx= NULL,
    spike.genome,
    spike.genome.idx= NULL,
    ref.gtf= NULL,
    fq.output.folder= "db/fq/PROSeq/",
    bam.output.folder= "db/bam/PROSeq/",
    alignment.stats.output.folder= "db/alignment_stats/PROSeq/",
    counts.output.folder= "db/counts/PROSeq/",
    counts.stats.output.folder= "db/stats/PROSeq/",
    bw.output.folder= "db/bw/PROSeq/",
    Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
    cores= 8,
    mem= 64,
    time= '2-00:00:00'
)
{
  # Retrieve ncRNA index ----
  if(!missing(ref.genome)) {
    ref.genome.ncRNA.idx <- switch(
      ref.genome,
      "mm10" = "/groups/stark/vloubiere/genomes/Mus_musculus/ncRNAs/bowtie_index/mm10_ncrna",
      NULL
    )
  }
  if(is.null(ref.genome.ncRNA.idx))
    stop("ref.genome.ncRNA.idx is missing, which is required to remove tRNAs/rRNAs.")

  # Trimming illumina adaptors ----
  cmd <- cmd_trimProseqAdaptors(fq1= fq1,
                                fq.output.folder= fq.output.folder)

  # * If several fq1 files provided, they will be merged during alignment ----
  fq1.trim <- cmd[file.type=="fq1.trim", path]

  # Align to rRNA and tRNA to remove them ----
  align.ncRNA.cmd <- cmd_alignBowtie(fq1= fq1.trim,
                                     output.prefix= output.prefix,
                                     genome.idx= ref.genome.ncRNA.idx,
                                     bam.output.folder= bam.output.folder,
                                     alignment.stats.output.folder = alignment.stats.output.folder,
                                     cores = cores)
  align.ncRNA.cmd[, file.type:= paste0(file.type, ".ncRNA")] # Make file types unique (see spike in below)
  cmd <- rbind(cmd, align.ncRNA.cmd, fill= TRUE)

  # Extract unaligned reads (valid reads) ----
  extract.cmd <- cmd_exractUnalignedReadsFromBam(bam = cmd[file.type=="bam.ncRNA", path],
                                                 fq.output.folder = fq.output.folder,
                                                 alignment.stats.output.folder = alignment.stats.output.folder,
                                                 cores = cores)
  extract.cmd[file.type=="fq1.unaligned", file.type:= "fq1.ncRNA.unaligned"]
  cmd <- rbind(cmd, extract.cmd, fill= TRUE)

  # Align to reference genome ----
  align.ref.cmd <- cmd_alignBowtie(fq1= cmd[file.type=="fq1.ncRNA.unaligned", path],
                                   output.prefix= output.prefix,
                                   genome= ref.genome,
                                   genome.idx= ref.genome.idx,
                                   bam.output.folder= bam.output.folder,
                                   alignment.stats.output.folder = alignment.stats.output.folder,
                                   cores = cores)
  align.ref.cmd[, file.type:= paste0(file.type, ".ref")] # Make file types unique (see ncRNA above)
  cmd <- rbind(cmd, align.ref.cmd, fill= TRUE)

  # Extract unaligned reads (spike-in) ----
  extract.cmd <- cmd_exractUnalignedReadsFromBam(bam = cmd[file.type=="bam.ref", path],
                                                 fq.output.folder = fq.output.folder,
                                                 alignment.stats.output.folder = alignment.stats.output.folder,
                                                 cores = cores)
  extract.cmd[file.type=="fq1.unaligned", file.type:= "fq1.ref.unaligned"]
  cmd <- rbind(cmd, extract.cmd, fill= TRUE)

  # Align spike-in reads (bowtie1) ----
  align.spike.cmd <- cmd_alignBowtie(fq1= cmd[file.type=="fq1.ref.unaligned", path],
                                     output.prefix= output.prefix,
                                     genome= spike.genome,
                                     genome.idx= spike.genome.idx,
                                     bam.output.folder= bam.output.folder,
                                     alignment.stats.output.folder = alignment.stats.output.folder,
                                     cores = cores)
  align.spike.cmd[, file.type:= paste0(file.type, ".spike")] # Make file types unique (see ref genome above)
  cmd <- rbind(cmd, align.spike.cmd, fill= TRUE)

  # Collapse and count total UMIs reference genome ----
  umi.ref.cmd <- cmd_umiCollapsingProseq(bam = cmd[file.type=="bam.ref", path],
                                         output.prefix = NULL, # From bam file
                                         counts.output.folder = counts.output.folder,
                                         stats.output.folder = counts.stats.output.folder)
  umi.ref.cmd[, file.type:= paste0(file.type, ".ref")] # Make file types unique (see spike in below)
  cmd <- rbind(cmd, umi.ref.cmd, fill= TRUE)

  # Collapse and count total UMIs spike-in genome ----
  umi.spike.cmd <- cmd_umiCollapsingProseq(bam = cmd[file.type=="bam.spike", path],
                                           output.prefix = NULL, # From bam file
                                           counts.output.folder = counts.output.folder,
                                           stats.output.folder = counts.stats.output.folder)
  umi.spike.cmd[, file.type:= paste0(file.type, ".spike")] # Make file types unique (see ref genome above below)
  cmd <- rbind(cmd, umi.spike.cmd, fill= TRUE)

  # Generate bw files ----
  bw.cmd <- cmd_umiToBigwigProseq(umi.counts = cmd[file.type=="umi.counts.ref", path],
                                  output.prefix = NULL, # From UMI counts file name
                                  bw.output.folder = bw.output.folder,
                                  Rpath = Rpath)
  cmd <- rbind(cmd, bw.cmd, fill= TRUE)

  # Return ----
  cmd$cores <- cores
  cmd$job.name <- paste0("PRO_", output.prefix)
  cmd$mem <- mem
  cmd$time <- time
  return(cmd)
}
