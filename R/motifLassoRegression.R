#' Perform LASSO Regression Using TF Motifs
#'
#' @description
#' This function performs LASSO regression to identify transcription factor (TF) motifs associated with a numeric response variable.
#' It splits the data into training and test sets, scales the data, performs cross-validation to select the best regularization parameter (`lambda`),
#' and fits a LASSO model. The function returns the model, test set predictions, and the most important predictors (motifs).
#'
#' @param response A numeric vector representing the variable to predict.
#' @param counts A matrix of motif counts, where rows correspond to observations and columns correspond to motifs.
#' @param names An optional character or factor vector of motif names. Defaults to the column names of `counts`.
#' @param CV A numeric value between 0 and 1 specifying the proportion of data to use as the test set. Default is `0.05`.
#' @param seed An integer specifying the random seed for reproducibility. Default is `123`.
#'
#' @return A list with the following components:
#' \item{test_set}{A `data.table` containing the observed and predicted values for the test set.}
#' \item{PCC}{The Pearson correlation coefficient between observed and predicted values on the test set.}
#' \item{R2}{The R-squared value (square of the Pearson correlation coefficient) for the test set.}
#' \item{best_predictors}{A `data.table` of the most important predictors (motifs) with non-zero coefficients, including their names and coefficients.}
#' \item{model}{The fitted LASSO model object.}
#'
#' @examples
#' # Download ATAC-Seq FC between PHD11 and control from the Nature paper
#' tmp <- tempfile(fileext = ".txt.gz")
#' download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222193/suppl/GSE222193%5FATAC%5FTransient%5Fph%2DKD%5FD11%5Fvs%5FControl%5Fno%5Fph%2DKD%5FDESeq2%5FFC.txt.gz",
#'               destfile = tmp)
#'
#' # Retrieve FC
#' FC <- fread(tmp)
#' FC <- cbind(importBed(FC$ID), FC[, !"ID"])
#'
#' # Resize
#' res <- resizeBed(bed = FC,
#'                  center = "center",
#'                  upstream = 250,
#'                  downstream = 250,
#'                  genome = "dm6")
#'
#' # Select JASPAR motifs
#' load("/groups/stark/vloubiere/motifs_db/vl_Dmel_motifs_DB_full.RData")
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds]
#' pwms <- vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_perc]
#' mot.cluster <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_cluster]
#'
#' # Compute motif counts
#' counts <- vl_motifCounts(bed= res,
#'                          genome= "dm6",
#'                          pwm_log_odds= sel)
#'
#' # Lasso regression
#' lasso <- motifLassoRegression(response= FC$log2FoldChange,
#'                               counts = counts,
#'                               names= mot.cluster)
#'
#' # Check PCC of the test set
#' vl_par()
#' plot(lasso$test_set)
#' addPcc(lasso$PCC)
#'
#' # Retrieve top predictors
#' preds <- lasso$best_predictors
#'
#' # Select top predictors
#' preds <- preds[, .SD[which.max(abs(s0))], name]
#' top <- preds[abs(s0)>0.03]
#' setorderv(top, "s0")
#'
#' # Plot
#' barplot(top$s0,
#'         horiz= TRUE,
#'         names.arg = top$name,
#'         xlab= "Lasso coefficient")
#'
#'
#' @export
motifLassoRegression <- function(response,
                                 counts,
                                 names = NULL,
                                 CV = 0.05,
                                 seed = 123) {
  # Checks
  if (is.data.table(counts))
    counts <- as.matrix(counts)
  if (!is.numeric(response))
    stop("Response must be numeric for regression.")
  if (length(response) != nrow(counts))
    stop("The number of rows in 'counts' must match the length of 'response'.")
  if (CV <= 0 || CV >= 1)
    stop("'CV' must be a value between 0 and 1.")
  if (is.null(names))
    names <- colnames(counts)
  if (!is.character(names) && !is.factor(names))
    stop("names should be a vector of characters or factors.")

  # Set random seed for reproducibility
  set.seed(seed)

  # Split data into training and test set
  test.set <- sample(c(TRUE, FALSE),
                     size = length(response),
                     replace = TRUE,
                     prob = c(CV, 1 - CV))
  testY <- response[test.set]
  testX <- counts[test.set, , drop = FALSE]
  print(paste(length(testY), "regions in test set."))
  trainY <- response[!test.set]
  trainX <- counts[!test.set, , drop = FALSE]
  print(paste(length(trainY), "regions in training set."))

  # Setting alpha = 1 implements lasso regression
  lambdas <- 10^seq(2, -3, by = -0.1)
  lasso_reg <- cv.glmnet(x = trainX,
                         y = trainY,
                         alpha = 1,
                         lambda = lambdas,
                         standardize = TRUE, # Let glmnet handle scaling
                         nfolds = 5,
                         family = "gaussian")
  # Best lambda
  lambda_best <- lasso_reg$lambda.min

  # Modelling
  model <- glmnet(x = trainX,
                  y = trainY,
                  alpha = 1,
                  lambda = lambda_best,
                  standardize = TRUE, # Let glmnet handle scaling
                  family = "gaussian")

  # Predicted values
  predict_test <- predict(model,
                          s = lambda_best,
                          newx = testX)
  predict_test <- as.vector(predict_test)

  # Test set
  test_set <- data.table(observed = testY,
                         predicted = predict_test)

  # Compute R-squared
  PCC <- cor(test_set$observed, test_set$predicted)
  rsq <- PCC^2

  # Retrieve best predictors
  preds <- as.matrix(model$beta)
  preds <- as.data.table(preds, keep.rownames = "motif")[s0 != 0]

  # Add names
  add.names <- data.table(motif = colnames(counts),
                          name = names)
  preds <- merge(preds, add.names)
  setorderv(preds, "s0")

  # Return
  return(list(test_set = test_set,
              PCC = PCC,
              R2 = rsq,
              best_predictors = preds,
              model = model))
}
