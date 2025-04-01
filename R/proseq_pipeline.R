#' Process PRO-Seq Sequencing Data
#'
#' @description
#' @export
proseqProcessing <- function(vbcFile,
                             layout,
                             eBC,
                             i7,
                             i5= 'none',
                             fq1,
                             fq2= NULL,
                             output.prefix,
                             ref.genome,
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
                             cores= 8)
{
  # Initiate command output ----
  cmd <- data.table(file.type= character(),
                    path= character(),
                    cmd= character())

  # Demultiplexing ----
  if(!missing(vbcFile)) {
    # Multiple bam/fq files can be provided for one sample
    .c <- data.table(vbcFile, i7, i5)
    message("The following samples will be demultiplexed:")
    message(paste(capture.output(print(.c)), collapse = "\n"))
    # Generate commands
    cmd <- .c[, {
      cmd_demultiplexVBCfile(vbcFile= vbcFile[1],
                             layout= layout,
                             i7= i7[1],
                             i5= i5[1],
                             output.prefix = output.prefix,
                             fq.output.folder= fq.output.folder,
                             start.seq = eBC[i],
                             trim.length = 10,
                             cores= cores)
    }, .(vbcFile, i7, i5)]
    cmd$vbcFile <- cmd$i7 <- cmd$i5 <- NULL
    # Set fq1 and fq2 parameters
    fq1 <- cmd[file.type=="fq1", path]
    if(layout=="PAIRED")
      fq2 <- cmd[file.type=="fq2", path]
  }

  # Trimming illumina adaptors ----
  # Multiple files can be provided for one sample
  for(i in seq(fq1)) {
    .c <- cmd_trimProseqAdaptors(fq1= fq1[i], # Second reads (fq2) are nots used
                                 fq.output.folder= fq.output.folder)
    cmd <- rbind(cmd, .c)
  }

  # * If several fq1 files provided, they will be merged at this step ----
  fq1.trim <- paste0(cmd[file.type=="fq1.trim", path], collapse = ",")

  # Alignment (bowtie1) ----
  align.ref.cmd <- cmd_alignBowtie(fq1= fq1.trim,
                                   output.prefix= output.prefix,
                                   genome= ref.genome,
                                   genome.idx= ref.genome.idx,
                                   mapq = NULL,
                                   bam.output.folder= bam.output.folder,
                                   alignment.stats.output.folder = alignment.stats.output.folder,
                                   cores = cores)
  align.ref.cmd[, file.type:= paste0(file.type, ".ref")] # Make file types unique (see spike in below)
  cmd <- rbind(cmd, align.ref.cmd)

  # Extract unaligned reads (spike-in) ----
  extract.cmd <- cmd_exractUnalignedReads(bam = bam.ref,
                                          fq.output.folder = fq.output.folder,
                                          alignment.stats.output.folder = alignment.stats.output.folder,
                                          cores = cores)
  cmd <- rbind(cmd, extract.cmd)

  # Align spike-in reads (bowtie1) ----
  align.spike.cmd <- cmd_alignBowtie(fq1= cmd[file.type=="fq1.unaligned", path],
                                     output.prefix= output.prefix,
                                     genome= spike.genome,
                                     genome.idx= spike.genome.idx,
                                     mapq = NULL,
                                     bam.output.folder= bam.output.folder,
                                     alignment.stats.output.folder = alignment.stats.output.folder,
                                     cores = cores)
  align.spike.cmd[, file.type:= paste0(file.type, ".spike")] # Make file types unique (see ref genome above)
  cmd <- rbind(cmd, align.spike.cmd)

  # UMI counts reference genome ----
  umi.ref.cmd <- cmd_umiToBigwigProseq(bam = cmd[file.type=="bam.ref", path],
                                       output.prefix = NULL, # From bam file
                                       counts.output.folder = counts.output.folder,
                                       stats.output.folder = counts.stats.output.folder)
  umi.ref.cmd[, file.type:= paste0(file.type, ".ref")] # Make file types unique (see spike in below)
  cmd <- rbind(cmd, umi.ref.cmd)

  # UMI counts spike-in genome ----
  umi.spike.cmd <- cmd_umiToBigwigProseq(bam = cmd[file.type=="bam.spike", path],
                                         output.prefix = NULL, # From bam file
                                         counts.output.folder = counts.output.folder,
                                         stats.output.folder = counts.stats.output.folder)
  umi.spike.cmd[, file.type:= paste0(file.type, ".spike")] # Make file types unique (see ref genome above below)
  cmd <- rbind(cmd, umi.spike.cmd)

  # Generate bw files ----
  bw.cmd <- cmd_umiToBigwigProseq(umi.counts = cmd[file.type=="umi.counts.ref"],
                                  output.prefix = NULL, # From UMI counts file name
                                  bw.output.folder = bw.output.folder,
                                  Rpath = Rpath)
  cmd <- rbind(cmd, bw.cmd)

  # Return ----
  return(cmd)
}
