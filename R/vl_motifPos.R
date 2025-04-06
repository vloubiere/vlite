#' Find motif positions
#'
#' Compute motif positions
#'
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified).
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start' and 'end' columns.
#' @param pwms_log_odds A PWMatrixList (in log2 odd ratio format) for which motif matches should be mapped. Overrides sel and motifDB arguments (see above).
#' @param genome Genome to be used for coordinates ("dm6, "dm3") and as background for counting motifs when bg= "genome".
#' @param bg Background used to find motifs. Possible values include "genome" and "even". Default= "genome"
#' @param p.cutoff p.value cutoff used for motif detection. For enrichment analyses based on presence/absence of a motif, high cutoff might perform better (1e-4 or 5e-5) while for regression analyses, lower cutoffs might be prefered (5e-4). Default= 5e-5 (stringent).
#' @param collapse.overlapping Should overlapping motifs be merged? If TRUE (default), motif instances that overlap more than 70 percent of their width are collapsed.
#' @param scratch Temp folder where temporary files will be saved. Useful to resume long jobs, see next argument
#' @param sub.folder Temp sub-folder where temporary files will be saved. Explicitly setting a sub.folder will allow the function to resume where it stoped at previous iteration. By default, sub.folder= tempdir().
#' @param overwrite If set to FALSE (default), motifs for which positions have already been computed and saved will be skipped.
#'
#' @examples
#' # Find position of 3 different motifs within two regions
#' pos <- vl_motifPos(vl_SUHW_top_peaks[1:2],
#'                     genome= "dm3",
#'                     pwm_log_odds= vl_Dmel_motifs_DB_full[c("cisbp__M2328", "flyfactorsurvey__suHw_FlyReg_FBgn0003567", "jaspar__MA0533.1"), pwms_log_odds, on= "motif_ID"])
#'
#' # Starting from sequence
#' pos <- vl_motifPos(sequence= "TGAGTTGTGTCTGAAATTGGGATTGCTGTTGCGACAATGCCTGTCTGACAGCATTGTCGATAAGAGCTTGAATCTGATTGGGGTCCATGGTAATATCTACCGTGGCACTATCTAACGGCCGACCTAATGCTTGGCCTACTTGCTCCTCCTCCCAGCTATCCTCGCTTTCGTATTCGACCTTAACCTTTCTGTAGTT#' ATGTGCCCAACTCATTGGTTGTTGGTTGGCACACCACAAATATACTGTTGCCGAGCACAATTGATCGGCTAAATGGTATGGCAAGAAAAGGTATGCAATATAATAATCTTTTATTGGGTATGCAACGAAAATTTGTTTCGTCAACGTATGCAATATTTTTTATTAAAAGAGGGTATGCAATGTATTTTATTAAAAACGGGTATGCAATATAATAATCTTTTATTGGG#' TATGCAACGAAAATTTGTTTCGTCAAAGTATGCAATATTTTTTATTAAAAGAGGGTATGCAATGTATTTTATTAAAAACGGGTATGCAATAAAAAATTATTTGGTTTCTCTAAAAAGTATGCAGCACTTATTTTTTGATAAGGTATGCAACAAAATTTTACTTTGCCGAAAATATGCAATGTTTTTGCGAATAAATTCAACGCACACTTATTACGTGGCCAGATACA#' CAACTTTTTTTTTTTTTTTTCACTCGTAAATTTCTTGATTGCGTCAAAGA",
#'                     genome= "dm3",
#'                     pwm_log_odds= vl_Dmel_motifs_DB_full[c("cisbp__M2328", "flyfactorsurvey__suHw_FlyReg_FBgn0003567", "jaspar__MA0533.1"), pwms_log_odds, on= "motif_ID"])
#'
#' # Motifs can also be mapped using a custom PWMatrixList, for example for promoter motifs
#' prom_db <- readRDS("/groups/stark/almeida/data/motifs/CP_motifs/CP_motifs_PWM.rds")
#' prom <- prom_db$Pwms_log_odds
#' for(i in seq(prom))
#'   prom[[i]]@profileMatrix <- pwmPercToLog(prom_db$Pwms_perc[[i]]@profileMatrix)
#'
#' # Select Ribo Protein gene promoters
#' proms <- rtracklayer::import("/groups/stark/annotations/dm3/dmel-all-filtered-r5.57_genes_and_transcripts_only.gff")
#' proms <- as.data.table(proms)
#' proms <- resizeBed(proms[grepl("ribosomal protein", fullname)], "start", 150, 50)
#' proms[, seqnames:= paste0("chr", seqnames)]
#'
#' prom_motifs <- vl_motifPos(proms,
#'                             pwm_log_odds= prom,
#'                             genome= "dm3")
#' # Many instances of TCT motif are found, as expected
#' prom_motifs$TC_17_Zabidi
#'
#' @return A list of positions of length = length(sequences)
#' @export
vl_motifPos <- function(sequences, ...) UseMethod("vl_motifPos")

#' @describeIn vl_motifPos Compute motifs position within
#' @export
vl_motifPos.data.table <- function(bed, genome, ...)
{
  sequences <- getBSsequence(bed, genome)
  vl_motifPos(sequences, genome= genome, ...)
}

#' @describeIn vl_motifPos Identify motif positions within sequences
#' @export
vl_motifPos <- function(sequences,
                        pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds],
                        genome,
                        bg= "genome",
                        p.cutoff= 5e-5,
                        collapse.overlapping= TRUE,
                        scratch= "/scratch/stark/vloubiere/",
                        sub.folder= tempdir(),
                        overwrite= FALSE)
{
  # Checks ----
  if(missing(genome))
    stop("genome is missing with no default")
  if(is.null(names(sequences)) || anyDuplicated(names(sequences)))
    stop("All sequences should have a unique name!")
  if(!"PWMatrixList" %in% class(pwm_log_odds))
    pwm_log_odds <- do.call(TFBSTools::PWMatrixList, pwm_log_odds)

  # Create tmp output folder ----
  tmp.folder <- paste0(scratch, "/", sub.folder, "/")
  print(paste0("Temp files will be stored in '", tmp.folder, "'"))
  dir.create(tmp.folder, showWarnings = FALSE, recursive = TRUE)

  # Metadata ----
  mot.names <- sapply(pwm_log_odds, TFBSTools::name)
  mot.names <- make.unique(mot.names)
  # File names cant handle '/'
  output <- paste0(tmp.folder, "/", gsub("/", "__", mot.names), ".rds")
  existing.files <- if(overwrite)
    rep(FALSE, length(output)) else
      file.exists(output)
  if(sum(existing.files))
  {
    print(paste0(sum(existing.files), "/", length(pwm_log_odds), " motif files already existed!"))
  }

  # Map motifs ----
  if(any(!existing.files))
  {
    parallel::mcmapply(function(mot, output.file)
    {
      res <- motifmatchr::matchMotifs(mot,
                                      sequences,
                                      genome= genome,
                                      p.cutoff= p.cutoff,
                                      bg= bg,
                                      out= "positions")[[1]]
      saveRDS(res, output.file)
    },
    mot= pwm_log_odds[!existing.files],
    output.file= output[!existing.files],
    mc.preschedule = TRUE,
    mc.cores = data.table::getDTthreads()-1)
  }

  # Processing ----
  pos <- data.table(motif= mot.names,
                    file= output)
  final <- pos[, {
    # Import
    .c <- readRDS(file)
    .c <- as.data.table(.c)
    # Add seqnames
    .c[, seqnames:= names(sequences)[group]]
    # If specified, collapsed motifs that overlap >70%, ignore.strand
    if(collapse.overlapping & nrow(.c))
    {
      min.gap <- -ceiling(mean(.c$width)*0.7)
      .c$idx <-  vl_collapseBed(.c,
                                min.gap = min.gap,
                                ignore.strand = TRUE,
                                return.idx.only = TRUE)
      .c <- .c[, {
        .(start= min(start),
          end= max(end),
          score= max(score, na.rm = TRUE))
      }, .(seqnames, idx)]
      .c[, width:= end-start+1]
    }
    # Simplify
    .c[, seqlvls:= seqnames]
    cols <- intersect(names(.c), c("start", "end", "strand", "width", "score"))
    .c <- .c[, .(mot.count= .N, ir= .(.SD)), seqlvls, .SDcols= cols]
    # Add missing levels
    all <- data.table(seqlvls= names(sequences))
    res <- merge(all, .c, by= "seqlvls", all.x= TRUE, sort= FALSE)
    # Return
    res
  }, motif]

  # Return ----
  return(final)
}
