#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Test if there are 9 args: if not, return an error ----
if (!length(args) %in% c(9, 11)) {
  stop(
    "Please specify:\n
       [required] 1/ The path to a .rds file containing SCT data for genes \n
       [required] 2/ The path to a .rds file containing SCT data for TFs (predictors) \n
       [required] 3/ The path to a .rds file containing SCT data clusters \n
       [required] 4/ Method. One of LASSO or ELASTIC \n
       [required] 5/ Number of iterations to perform \n
       [required] 6/ The fraction of cells to use per iteration \n
       [required] 7/ The fraction of times a TF should be recovered to be a candidate \n
       [required] 8/ Output prefix \n
       [required] 9/ Output folder \n
       [optional] 10/ Starting line index for the genes matrix \n
       [optional] 11/ Ending line index for the genes matrix \n"
  )
}

# Load libraries ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Matrix, warn.conflicts = FALSE))
suppressMessages(library(glmnet, warn.conflicts = FALSE))
suppressMessages(library(parallel, warn.conflicts = FALSE))

# args <- c("/scratch-cbe/users/vincent.loubiere/candidate_regulons/82bdf52dd5f421b69e97350f873121e2_genes_norm_counts.rds",
#           "/scratch-cbe/users/vincent.loubiere/candidate_regulons/ff1e664a1fe0a395a28705f6760b6b45_TFs_norm_counts.rds",
#           "/scratch-cbe/users/vincent.loubiere/candidate_regulons/d11c7fada862b5ffea2d6a64f8017627_cell_labels.rds",
#           "ELASTIC",
#           "5",
#           "0.3",
#           "0.6",
#           "genes_1_20",
#           "db/candidate_regulons/",
#           "1",
#           "20")

# Parse arguments ----
SCT_data_genes <- args[1]
SCT_data_TFs <- args[2]
cell_labels <- args[3]
method <- args[4]
N_iter <- as.numeric(args[5]) # Number of iterations. > 50
frac_cells <- as.numeric(args[6]) # Number of cells per iteration. Generally .3
stability <- as.numeric(args[7]) # Fraction required to be considered as candidate. Generally .6
output_prefix <- args[8]
output_folder <- args[9]
if(length(args==11)) {
  starting.line <- as.integer(args[10])
  ending.line <- as.integer(args[11])
}
stopifnot(grepl(".rds$", SCT_data_genes))
stopifnot(grepl(".rds$", SCT_data_TFs))
stopifnot(grepl(".rds$", cell_labels))
stopifnot(method %in% c("LASSO", "ELASTIC"))

# Import genes and TF expression matrices ----
genes <- readRDS(SCT_data_genes)
if(length(args==11))
  genes <- genes[starting.line:ending.line,]
TFs <- readRDS(SCT_data_TFs)
stopifnot(class(genes)=="dgCMatrix")
stopifnot(class(TFs)=="dgCMatrix")
stopifnot(ncol(genes)==ncol(TFs))

# Import cell labels ----
labels <- readRDS(cell_labels)
stopifnot(is.data.table(labels))
stopifnot(all(c("cluster", "cdition") %in% names(labels)))
clusters <- labels$cluster
cditions <- factor(labels$cdition)  # ensure factor with desired baseline level
cdition_mat <- model.matrix(~ cditions)[, -1, drop = FALSE]  # C-1 dummies, no intercept
# cdition_mat is n_cells x (C-1). The first condition will be used as base effect (W18)

# Compute clusters probabilities ----
probs <- (table(clusters)/length(clusters))[clusters]

# Loop over genes ----
final <- parallel::mclapply(
  rownames(genes),
  function(x) {
    # Loop over iterations ----
    t0 <- Sys.time()
    res <- lapply(seq(N_iter), function(i) {
      set.seed(i)
      train_sel <- seq(ncol(genes)) %in% sample(ncol(genes), round(frac_cells*ncol(genes)), prob = probs)

      # Base TF matrix excluding the target TF itself (if present)
      TF_keep <- rownames(TFs) != x
      X_TF_train <- t(TFs[TF_keep, train_sel, drop = FALSE])    # n_train x p
      X_TF_full  <- t(TFs[TF_keep, , drop = FALSE])             # n_all x p

      # Condition dummies (align rows to cells)
      X_cond_train <- cdition_mat[train_sel, , drop = FALSE]       # n_train x (C-1)
      X_cond_full  <- cdition_mat                                   # n_all x (C-1)

      # Interactions: column-wise multiply each TF column by each condition column
      make_interactions <- function(X_tfs, X_cond) {
        # Kronecker-style column interactions: for each cond col, multiply all TF cols
        # Result: n x (p * (C-1))
        do.call(cbind, lapply(seq_len(ncol(X_cond)), function(j) X_tfs * X_cond[, j]))
      }
      X_int_train <- make_interactions(X_TF_train, X_cond_train)
      X_int_full  <- make_interactions(X_TF_full,  X_cond_full)

      # Combine: [TF main | condition | TF×condition]
      X_train <- cbind(X_TF_train, X_cond_train, X_int_train)
      X_full  <- cbind(X_TF_full,  X_cond_full,  X_int_full)

      # Response
      y_train <- genes[x, train_sel]
      y_full  <- genes[x, ]

      # Modelling using elastic net or lasso
      fit <- if(method=="ELASTIC") {
        # Grid of alphas (try a few; 0.7–0.9 often works well)
        alpha_grid <- c(0.3, 0.5, 0.7, 0.9)

        best <- list(cvm = Inf, alpha = NA, lambda = NA, cvfit = NULL)
        lambdas <- 10^seq(2, -3, by = -0.1)

        for (a in alpha_grid) {
          cvfit <- glmnet::cv.glmnet(
            x = X_train, y = y_train,
            alpha = a,
            lambda = lambdas,
            standardize = TRUE,
            nfolds = 5,
            family = "gaussian"
          )
          if (min(cvfit$cvm) < best$cvm) {
            best$cvm <- min(cvfit$cvm)
            best$alpha <- a
            best$lambda <- cvfit$lambda.min
            best$cvfit <- cvfit
          }
        }

        # Fit model
        glmnet::glmnet(
          x = X_full,
          y = y_full,
          alpha = best$alpha,
          lambda = best$lambda,
          standardize = TRUE,
          family = "gaussian"
        )
      } else if(method=="LASSO") {
        # LASSO with CV to pick lambda
        lambdas <- 10^seq(2, -3, by = -0.1)
        cvfit <- glmnet::cv.glmnet(
          x = X_train,
          y = y_train,
          alpha = 1,
          lambda = lambdas,
          standardize = TRUE,
          nfolds = 5,
          family = "gaussian"
        )
        lambda_best <- cvfit$lambda.min

        # Fit model
        glmnet::glmnet(
          x = X_full,
          y = y_full,
          alpha = 1,
          lambda = lambda_best,
          standardize = TRUE,
          family = "gaussian"
        )
      }

      # Name coefficients: TF, condition, interactions
      TF_names <- rownames(TFs)[TF_keep]
      cond_names <- colnames(X_cond_full)  # like conditionB, conditionC...
      int_names <- as.vector(outer(TF_names, cond_names, paste, sep="__x__"))
      colnames_all <- c(TF_names, cond_names, int_names)
      coefs <- as.matrix(fit$beta)
      rownames(coefs) <- colnames_all

      # Return as data.table with columns: term, coef
      return(data.table::data.table(term = rownames(coefs), s0 = as.numeric(coefs[, 1])))
    })
    t1 <- Sys.time()
    # print(t1-t0)

    # Cbind ----
    mat <- do.call(cbind, lapply(res, function(x) x$s0))

    # Aggregate results ----
    agg <- res[[1]][, .(term)]
    agg[, c("TF", "cdition"):= tstrsplit(term, "__x__cditions")]
    agg[is.na(cdition), cdition:= "base"]

    # Remove conditions coeff (can be ignored to infer regulons) ----
    keep <- !grepl("^cditions", agg$term)
    mat <- mat[(keep),]
    agg <- agg[(keep)]

    # # Add base effect (W18) to all condition-specific coeffs (commented out for now, keep raw coeffs) ----
    cdition.mat <- mat
    # base.idx <- agg$cdition == "base"
    # for(cdition in setdiff(unique(agg$cdition), "base")) {
    #   cdition.idx <- agg$cdition == cdition
    #   cdition.mat[cdition.idx,] <- mat[base.idx,]+mat[cdition.idx,]
    # }

    # Compute metrix ----
    agg[, pos.frac:= apply(cdition.mat, 1, function(x) sum(x>0)/length(x))]
    agg[, neg.frac:= apply(cdition.mat, 1, function(x) sum(x<0)/length(x))]
    agg[, median.abs:= apply(cdition.mat, 1, function(x) median(abs(x[x!=0])))]
    agg[, median.signed:= apply(cdition.mat, 1, function(x) median(x[x!=0]))]
    # Candidate regulosn
    agg[, act.cand:= pos.frac>=stability]
    agg[, rep.cand:= neg.frac>=stability]
    agg$term <- NULL

    # Return ----
    return(agg)
  },
  mc.preschedule = TRUE,
  mc.cores = data.table::getDTthreads()-1
)

# Save ----
names(final) <- rownames(genes)
final <- rbindlist(final, idcol = "gene_symbol")
saveRDS(final,
        file.path(output_folder, paste0(output_prefix, "_candidate_regulons_", method, ".rds")))
