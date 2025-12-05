#' Compute candidate regulons for all genes
#'
#' Wrapper around the cmd_infer_candidate_regulons command.
#'
#' @param seurat.object
#' @param N.variable.genes
#' @param min.act.cells
#' @param TF.gene.names
#' @param N.genes
#' @param cores
#' @param mem
#' @param output.prefix
#' @param output.folder
#' @param tmp.folder
#' @param overwrite
#' @param N.iterations
#' @param cell.frac
#' @param stability.cutoff
#' @param method
#'
#' @return
#' @export
#'
#' @examples
SCENIClite_infer_regulons <- function(seurat.object,
                                      N.variable.genes= NULL,
                                      min.act.cells= 0.025,
                                      TF.gene.names,
                                      N.genes= cores,
                                      N.iterations= 100L,
                                      cell.frac= .3,
                                      stability.cutoff= .6,
                                      method= "ELASTIC",
                                      output.prefix,
                                      output.folder= "db/candidate_regulons/",
                                      cores= 20,
                                      mem= cores*2,
                                      tmp.folder= "/scratch-cbe/users/vincent.loubiere/candidate_regulons/",
                                      overwrite= F)
{
  # Checks
  if(is.character(seurat.object) && grepl(".rds$", seurat.object))
    seurat.object <- readRDS(seurat.object)
  stopifnot(class(seurat.object)[1] == "Seurat")
  stopifnot("sample_id" %in% names(dat@meta.data)) # sample-aware model
  stopifnot(length(dat@assays$SCT@var.features)>0 || is.null(N.variable.genes)) # Variable genes need to be defined
  stopifnot(is.character(TF.gene.names))
  stopifnot(is.numeric(N.genes) || N.genes<100) # This would get very long
  stopifnot(min.act.cells<0.1) # This would get very long

  # Final output file ----
  output.file <- file.path(output.folder, paste0(output.prefix, "_candidate_regulons.rds"))

  # Cached input files ----
  cell.labels.file <- vl_cache_file(
    input.list = list(rownames(dat), colnames(dat)),
    prefix = "cell_labels",
    tmp.dir = tmp.folder
  )
  genes.file <- vl_cache_file(
    input.list = list(cell.labels.file, N.variable.genes, min.act.cells),
    prefix = "genes_norm_counts",
    tmp.dir = tmp.folder
  )
  TFs.file <- vl_cache_file(
    input.list = list(genes.file, TF.gene.names),
    prefix = "TFs_norm_counts",
    tmp.dir = tmp.folder
  )

  # If the command is to be executed ----
  if(overwrite | !file.exists(output.file))
  {
    # Create temporary dir ----
    dir.create(tmp.folder, recursive = T, showWarnings = F)

    # Save genes matrix ----
    if(!file.exists(genes.file) | overwrite) {
      # Extract variable genes
      if(!is.null(N.variable.genes))
        dat <- Seurat::FindVariableFeatures(dat, nfeatures= N.variable.genes)
      SCT <- dat[["SCT"]]$data[VariableFeatures(dat),]
      # Activity filter
      Ncells <- round(ncol(SCT)*min.act.cells)
      sel <- Matrix::rowSums(SCT>0) >= Ncells
      print(paste0(sum(sel), "/", nrow(SCT), " variable genes are active in at least ", Ncells, " cells and will be used."))
      SCT <- SCT[sel,]
      # Save
      saveRDS(SCT, genes.file)
    } else
      SCT <- readRDS(genes.file)

    # Save TF predictors matrix ----
    if(!file.exists(TFs.file) | overwrite) {
      # Retrieve TF genes
      TF.gene.names <- unique(TF.gene.names)
      sel <- intersect(TF.gene.names, rownames(SCT))
      print(paste0(length(sel), "/", length(TF.gene.names), " TF genes were found and will be used as predictors."))
      TFs <- SCT[sel,]
      # Save
      saveRDS(TFs, TFs.file)
    }

    # Save labels and sample IDs ----
    if(!file.exists(cell.labels.file) | overwrite) {
      labels <- data.table(
        cluster= Seurat::Idents(dat),
        cdition= dat$sample_id
      )
      # Save
      saveRDS(labels, cell.labels.file)
    }

    # Split jobs ----
    cmd <- data.table(starting.line= seq(1, nrow(SCT), N.genes))
    cmd[, ending.line:= starting.line+N.genes-1]
    cmd[.N, ending.line:= nrow(SCT)]

    # Generate commands ----
    cmd <- cmd[, {
      cmd_infer_candidate_regulons(
        genes.file = genes.file,
        TFs.file = TFs.file,
        cell.labels.file = cell.labels.file,
        N.iterations = N.iterations,
        cell.frac = cell.frac,
        stability.cutoff = stability.cutoff,
        method = method,
        output.folder = paste0(tmp.folder, "/tmp/"),
        starting.line= starting.line,
        ending.line= ending.line
      )
    }, .(starting.line, ending.line)]

    # Submit split jobs to the server ----
    if(overwrite | any(!file.exists(cmd$path))) {
      print(paste(sum(!file.exists(cmd$path)), "jobs with", N.genes, " genes each will now be submitted:"))

      # Submit commands
      cmd[, {
        vl_submit(
          .SD,
          execute = TRUE,
          overwrite = overwrite,
          logs= paste0(tmp.folder, "/logs/")
        )
      }, starting.line]

    } else {
      # Aggregate results
      cand <- rbindlist(lapply(cmd$path, readRDS))
      # Save
      dir.create(output.folder, showWarnings = F, recursive = T)
      saveRDS(cand, output.file)
    }
  } else {
    message(paste0("Reading output file -> ", output.file))
    readRDS(output.file)
  }
}
