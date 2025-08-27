#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Test if there are 7 args: if not, return an error ----
if (length(args) != 7) {
  stop(
    "Please specify:\n
       [required] 1/ The path to a .rds file containing SCT data for genes \n
       [required] 2/ The path to a .rds file containing SCT data for TFs (predictors) \n
       [required] 3/ The path to a .rds file containing SCT data clusters \n
       [required] 4/ Number of iterations to perform \n
       [required] 5/ The fraction of cells to use per iteration \n
       [required] 6/ The fraction of times a TF should be recovered to be a candidate \n
       [required] 7/ Output folder \n"
  )
}

# Load libraries ----
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(Matrix, warn.conflicts = FALSE))
suppressMessages(library(glmnet, warn.conflicts = FALSE))
suppressMessages(library(parallel, warn.conflicts = FALSE))

# Parse arguments ----
SCT_data_genes <- args[1]
SCT_data_TFs <- args[2]
SCT_data_clusters <- args[3]
N_iter <- as.numeric(args[4]) # Number of iterations. > 50
frac_cells <- as.numeric(args[5]) # Number of cells per iteration. Generally .3
stability <- as.numeric(args[6]) # Fraction required to be considered as candidate. Generally .6
output_folder <- args[7]
stopifnot(grepl(".rds$", SCT_data_genes))
stopifnot(grepl(".rds$", SCT_data_TFs))
stopifnot(grepl(".rds$", SCT_data_clusters))

# Parse arguments ----
genes <- readRDS(SCT_data_genes)
TFs <- readRDS(SCT_data_TFs)
stopifnot(class(genes)=="dgCMatrix")
stopifnot(class(TFs)=="dgCMatrix")
stopifnot(ncol(genes)==ncol(TFs))
clusters <- readRDS(SCT_data_clusters)
stopifnot(is.data.table(clusters))
stopifnot(ncol(clusters)==1)
clusters <- clusters[[1]]

# Compute clusters probabilities ----
probs <- (table(clusters)/length(clusters))[clusters]

# Loop over genes ----
parallel::mclapply(
  rownames(genes),
  function(x) {
    # Loop over iterations ----
    res <- lapply(seq(N_iter), function(i) {
      set.seed(i)
      train_sel <- seq(ncol(genes)) %in% sample(ncol(genes), round(frac_cells*ncol(genes)), prob = probs)

      # Setting alpha = 1 implements lasso regression
      lambdas <- 10^seq(2, -3, by = -0.1)
      lasso_reg <- glmnet::cv.glmnet(x = t(TFs[rownames(TFs) != x, train_sel]),
                                     y = genes[x, train_sel],
                                     alpha = 1,
                                     lambda = lambdas,
                                     standardize = TRUE, # Let glmnet handle scaling
                                     nfolds = 5,
                                     family = "gaussian")

      # Best lambda
      lambda_best <- lasso_reg$lambda.min

      # Modelling
      model <- glmnet::glmnet(x = t(TFs[rownames(TFs) != x, ]),
                              y = genes[x, ],
                              alpha = 1,
                              lambda = lambda_best,
                              standardize = TRUE, # Let glmnet handle scaling
                              family = "gaussian")
      # Return
      return(as.data.table(as.matrix(model$beta), keep.rownames = "symbol"))
    })

    # Aggregate results ----
    agg <- cbind(res[[1]], do.call(cbind, lapply(res[-1], function(x) x[, .(s0)])))
    act_cand <- rowSums(agg[,-1]>0) >= stability*N_iter
    rep_cand <- rowSums(agg[,-1]<0) >= stability*N_iter
    median_abs <- apply(agg[,-1], 1, function(x) median(abs(x[x!=0])))
    median_signed <- apply(agg[,-1], 1, function(x) median(x[x!=0]))
    agg[, act_cand:= act_cand]
    agg[, rep_cand:= rep_cand]
    agg[, median_abs:= median_abs]
    agg[, median_signed:= median_signed]
    names(agg) <- make.unique(names(agg))
    setcolorder(agg,
                c("symbol", "act_cand", "rep_cand", "median_abs", "median_signed"))

    # Save ----
    saveRDS(agg,
            file.path(output_folder, paste0(x, "_lasso.rds")))
  },
  mc.preschedule = TRUE,
  mc.cores = data.table::getDTthreads()-1
)
