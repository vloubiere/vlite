#' Find Motif Positions in Sequences or Genomic Regions
#'
#' @description
#' Identifies the positions of motifs in a set of sequences or genomic regions using a PWMatrixList.
#'
#' @param sequences A named character vector of sequences to analyze. This argument takes precedence over the bed argument.
#' @param bed Genomic ranges in a format compatible with ?importBed, from which genomic sequences will be retrieved when sequences is set to NULL.
#' @param pwm_log_odds A PWMatrixList (in log2 odds ratio format) containing motifs to map. For example, see "/groups/stark/vloubiere/motifs_db/".
#' @param genome The genome to use as background when bg = "genome" and/or to retrieve sequences (when bed is specified). This argument is required.
#' @param bg The background model for motif detection. Options are "genome", "subject" (inferred from input sequences) or "even" (0.25, 0.25, 0.25, 0.25). Default= "genome".
#' @param p.cutoff The p-value cutoff for motif detection. Default= 5e-5.
#' @param pos.strand If set to TRUE, only motifs on the positive strand are considered (default= FALSE).
#' @param collapse.overlapping Logical. If TRUE, motifs with 70 percent overlap ore more will be collapsed, irrespective of their strand. Default= FALSE.
#' @param scratch The path to the scratch directory for temporary files. Default= "/scratch-cbe/users/vincent.loubiere/motifs/".
#' @param cleanup.cache If set to TRUE, overwrites cached intermediate results. Default= FALSE.
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
#' pos1 <- vl_motifPos(bed= suhw,
#'                     genome= "dm3",
#'                     pwm_log_odds= sel,
#'                     overwrite= FALSE)
#'
#' # Starting from sequence
#' seq <- getBSsequence(suhw, "dm3")
#' pos2 <- vl_motifPos(sequences= seq,
#'                     genome= "dm3",
#'                     pwm_log_odds= sel,
#'                     overwrite= FALSE)
#'
#' # Should give similar results
#' identical(pos1, pos2)
#'
#' @export
vl_motifPos <- function(sequences,
                        bed,
                        pwm_log_odds,
                        genome,
                        bg= "genome",
                        p.cutoff= 5e-5,
                        pos.strand= FALSE,
                        collapse.overlapping= FALSE,
                        scratch= "/scratch-cbe/users/vincent.loubiere/motifs/",
                        cleanup.cache= FALSE)
{
  # Checks ----
  if(!missing(sequences) && !missing(bed))
    warning("sequences are provided -> input bed will not be used.")
  if(missing(sequences) && missing(bed))
    stop("sequences of bed regions should be specified.")
  if(missing(genome) && (missing(sequences) | bg=="genome"))
    stop("genome is missing with no default.")
  if(!"PWMatrixList" %in% class(pwm_log_odds))
    pwm_log_odds <- do.call(TFBSTools::PWMatrixList, pwm_log_odds)

  # Get sequences ----
  if(missing(sequences)) {
    bed <- importBed(bed)
    sequences <- getBSsequence(bed, genome)
  }

  # Make sure sequence names are unique ----
  if(is.null(names(sequences)) || anyDuplicated(names(sequences)))
    stop("All sequences should have a unique name!")

  # Create tmp output folder ----
  params <- list(sequences,
                 ifelse(bg=="genome", genome, bg),
                 bg,
                 p.cutoff)
  sub.folder <- digest::digest(params)
  tmp.folder <- paste0(scratch, "/", sub.folder, "/")
  print(paste0("Temp files will be stored in '", tmp.folder, "'"))
  dir.create(tmp.folder,
             showWarnings = FALSE,
             recursive = TRUE)

  # Create final file cache name ----
  params <- list(sub.folder,
                 pwm_log_odds,
                 pos.strand,
                 collapse.overlapping)
  final.file <- paste0(tmp.folder, digest::digest(params), ".rds")

  # Main function ----
  if(!file.exists(final.file) | cleanup.cache) {

    # Retrieve unique motif names ----
    mot.names <- sapply(pwm_log_odds, TFBSTools::name)
    mot.names <- make.unique(mot.names)

    # Generate output file names ----
    output.files <- paste0(tmp.folder, "/", gsub("/", "__", mot.names), ".rds")

    # Check for existing files ----
    missing.files <- !file.exists(output.files) | cleanup.cache

    # Map motifs ----
    if(any(missing.files))
    {
      print(paste0(sum(!missing.files), "/", length(pwm_log_odds), " motif files already existed!"))

      parallel::mcmapply(function(mot, output.file)
      {
        res <- motifmatchr::matchMotifs(pwms = mot,
                                        subject = sequences,
                                        genome = genome,
                                        p.cutoff = p.cutoff,
                                        bg = bg,
                                        out = "positions")[[1]]
        saveRDS(res, output.file)
      },
      mot= pwm_log_odds[missing.files],
      output.file= output.files[missing.files],
      mc.preschedule = TRUE,
      mc.cores = max(c(1, data.table::getDTthreads()-1)))
    }
    print("All motif positions computed ;)")

    # Processing ----
    pos <- data.table(motif= mot.names,
                      file= output.files)
    final <- pos[, {
      # Import
      .c <- readRDS(file)
      .c <- as.data.table(.c)
      # Select positive strand motifs
      if(pos.strand)
        .c <- .c[strand=="+"]
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

    # Save cache file
    saveRDS(final,
            final.file)
  } else
    final <- readRDS(final.file)

  # Return ----
  return(final)
}
