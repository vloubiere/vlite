#' Find Motif Positions in Sequences or Genomic Regions
#'
#' @description
#' Identifies the positions of motifs in a set of sequences or genomic regions using a PWMatrixList.
#'
#' @param sequences A named character vector of sequences to analyze. This argument takes precedence over the bed argument.
#' @param bed Genomic ranges in a format compatible with ?importBed(). Used to retrieve sequences if sequences is not provided.
#' @param pwm_log_odds A PWMatrixList (in log2 odds ratio format) containing motifs to map. For example, see "/groups/stark/vloubiere/motifs_db/".
#' @param genome The genome to use as background when bg = "genome" and/or to retrieve sequences (when bed is specified). This argument is required.
#' @param bg The background model for motif detection. Options are "genome" or "even". Default= "genome".
#' @param p.cutoff The p-value cutoff for motif detection. Default= 5e-5.
#' @param collapse.overlapping Logical. If TRUE, overlapping motifs (greater than 70 percent overlap) are collapsed into a single region. Default= TRUE.
#' @param scratch The path to the scratch directory for temporary files. Default= "/scratch-cbe/users/vincent.loubiere/motifs/".
#' @param sub.folder A subfolder within the scratch directory for temporary files. If set to NULL (default), a unique id will be generated using digest::digest(list(function.parameters)).
#' @param overwrite If set to TRUE, overwrites cached intermediate results. Default= FALSE.
#'
#' @return A data.table containing motif positions for each sequence, with the following columns:
#' - `motif`: The name of the motif.
#' - `seqlvls`: The sequence names.
#' - `mot.count`: The number of motifs detected.
#' - `ir`: A nested `data.table` with columns for motif positions (`start`, `end`, `strand`, `width`, `score`).
#'
#' @examples
#' # Dm3 coordinates of the two strongest SU(HW) peaks
#' suhw <- importBed(c("chr2R:942901-943600", "chrX:1417201-1418000"))
#'
#' # Load SU(HW) motifs
#' load("/groups/stark/vloubiere/motifs_db/vl_Dmel_motifs_DB_full.RData")
#' sel <- vl_Dmel_motifs_DB_full[c("cisbp__M2328", "flyfactorsurvey__suHw_FlyReg_FBgn0003567", "jaspar__MA0533.1"), pwms_log_odds, on= "motif_ID"]
#'
#' # Find position of 3 different motifs within two regions
#' pos1 <- vl_motifPos.data.table(bed= suhw,
#'                                genome= "dm3",
#'                                pwm_log_odds= sel,
#'                                overwrite= FALSE)
#'
#' # Starting from sequence
#' seq <- getBSsequence(suhw, "dm3")
#' pos2 <- vl_motifPos(sequence= seq,
#'                     genome= "dm3",
#'                     pwm_log_odds= sel,
#'                     overwrite= FALSE)
#'
#' # Should give similar results
#' identical(pos1, pos2)
#'
#'
#' @export
vl_motifPos <- function(sequences = NULL, bed = NULL, ...) {
  if (!is.null(bed)) {
    if (!is.data.table(bed)) {
      bed <- importBed(bed)
    }
    return(vl_motifPos.data.table(bed = bed, ...))
  }
  vl_motifPos.default(sequences = sequences, ...)
}

#' @describeIn vl_motifPos Compute motifs position within
#' @export
vl_motifPos.data.table <- function(bed, genome, ...)
{
  sequences <- getBSsequence(bed, genome)
  vl_motifPos.default(sequences= sequences,
                      genome= genome,
                      ...)
}

#' @describeIn vl_motifPos Identify motif positions within sequences
#' @export
vl_motifPos.default <- function(sequences,
                                pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds],
                                genome,
                                bg= "genome",
                                p.cutoff= 5e-5,
                                collapse.overlapping= TRUE,
                                scratch= "/scratch-cbe/users/vincent.loubiere/motifs/",
                                sub.folder= NULL,
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
  if(is.null(sub.folder)) {
    params <- list(sequences,
                   pwm_log_odds,
                   genome,
                   bg,
                   p.cutoff)
    sub.folder <- digest::digest(params)
  }
  tmp.folder <- paste0(scratch, "/", sub.folder, "/")
  print(paste0("Temp files will be stored in '", tmp.folder, "'"))
  dir.create(tmp.folder, showWarnings = FALSE, recursive = TRUE)

  # Retrieve unique motif names ----
  mot.names <- sapply(pwm_log_odds, TFBSTools::name)
  mot.names <- make.unique(mot.names)

  # Generate output file names ----
  output <- paste0(tmp.folder, "/", gsub("/", "__", mot.names), ".rds")

  # Check for existing files ----
  existing.files <- if(overwrite) {
    rep(FALSE, length(output))
  } else {
    file.exists(output)
  }
  if(sum(existing.files)) {
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
      .c$idx <-  collapseBed(bed = .c,
                             max.gap = min.gap,
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
