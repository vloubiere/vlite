#' Process ORtag Sequencing Data
#'
#' @param vbcFile
#' @param layout
#' @param i7
#' @param i5
#' @param fq1
#' @param output.prefix
#' @param genome
#' @param genome.idx
#' @param gtf
#' @param fq.output.folder
#' @param bam.output.folder
#' @param alignment.stats.output.folder
#' @param bed.output.folder
#' @param counts.output.folder
#' @param Rpath
#' @param cores
#'
#' @description
#' @export
orftagProcessing <- function(vbcFile,
                             layout,
                             i7,
                             i5= 'none',
                             fq1,
                             output.prefix,
                             genome,
                             genome.idx= NULL,
                             gtf= NULL,
                             fq.output.folder= "db/fq/ORFtag/",
                             bam.output.folder= "db/bam/ORFtag/",
                             alignment.stats.output.folder= "db/alignment_stats/ORFtag/",
                             bed.output.folder= "db/bed/ORFtag/",
                             counts.output.folder= "db/counts/ORFtag/",
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
                             cores= cores)
    }, .(vbcFile, i7, i5)]
    cmd$vbcFile <- cmd$i7 <- cmd$i5 <- NULL
    # Set fq1 parameter
    fq1 <- cmd[file.type=="fq1", path]
  }

  # Trimming illumina adaptors ----
  # Multiple files can be provided for one sample
  for(i in seq(fq1)) {
    .c <- cmd_trimIlluminaAdaptors(fq1= fq1[i],
                                   fq2= NULL, # Second reads are not used
                                   layout= "SINGLE",
                                   fq.output.folder= fq.output.folder)
    cmd <- rbind(cmd, .c)
  }

  # * If several fq1 files provided, they will be merged at this step ----
  fq1.trim <- paste0(cmd[file.type=="fq1.trim", path], collapse = ",")

  # Alignment ----
  align.cmd <- cmd_alignBowtie2(fq1= fq1.trim,
                                fq2= NULL, # Second reads are not used
                                layout= "SINGLE",
                                output.prefix= output.prefix,
                                genome= genome,
                                genome.idx= genome.idx,
                                mapq= 30,
                                bam.output.folder= bam.output.folder,
                                alignment.stats.output.folder= alignment.stats.output.folder,
                                cores= cores)
  cmd <- rbind(cmd, align.cmd)

  # Collapse bam file (unique insertions) ----
  collapse.cmd <- cmd_collapseBam(bam = cmd[file.type=="bam", path],
                                  output.prefix = NULL, # From bam file
                                  collapsed.bam.output.folder = bam.output.folder,
                                  collapsed.stats.output.folder = alignment.stats.output.folder,
                                  cores = cores)
  cmd <- rbind(cmd, collapse.cmd)

  # Assign insertions to closest downstream exon ----
  cmd <- cmd_assignInsertions(bam = cmd[file.type=="collapsed.bam", path],
                              output.prefix = NULL, # From bam file
                              genome = genome,
                              gtf = gtf,
                              bed.output.folder = bed.output.folder,
                              counts.output.folder = counts.output.folder,
                              Rpath = Rpath)
  cmd <- rbind(cmd, collapse.cmd)

  # Return ----
  return(cmd)
}
