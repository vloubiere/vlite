#' Perform Regularized Regression (LASSO / Elastic Net) Using TF Motifs
#'
#' @description
#' Fits a regularized linear model (glmnet) to predict a numeric response from motif-count features.
#' Supports LASSO (alpha = 1) and Elastic Net (0 < alpha < 1). Repeats random train/test splits,
#' does inner CV on the training set to pick lambda, fits the model, and summarizes stability of selected motifs.
#'
#' @param response Numeric vector to predict (length = nrow(counts)).
#' @param counts Matrix of motif counts (n_regions x n_motifs).
#' @param names Optional motif display names (length = ncol(counts)). Defaults= colnames(counts).
#' @param CV Test set proportion in (0,1). Default= 0.1.
#' @param seed Integer seed. Default 123.
#' @param repeats Integer number of repeats. Default= 10.
#' @param alpha Numeric in [0,1]. alpha=1 is LASSO; 0<alpha<1 is Elastic Net; alpha=0 is ridge. Default= 0.5.
#' @param use_1se Logical; if TRUE uses lambda.1se instead of lambda.min. Default= TRUE (more stable).
#' @param nfolds Integer folds for inner CV. Default= 5.
#'
#' @return List with:
#' \item{stability}{data.table summarizing selection frequency and coefficient stability.}
#' \item{performance}{data.table with PCC/R2/lambda per repeat.}
#' \item{test_set}{data.table with observed/predicted across repeats.}
#' \item{best_predictors}{data.table of non-zero coefficients per repeat.}
#' \item{models}{list of fitted glmnet models (one per repeat).}
#' \item{alpha}{The alpha value that was used.}
#' \item{use_1se}{use_1se value.}
#'
#' @export
motifRegularizedRegression <- function(response,
                                       counts,
                                       names = NULL,
                                       CV = 0.1,
                                       seed = 123,
                                       repeats = 10,
                                       alpha = 0.5,
                                       use_1se = TRUE,
                                       nfolds = 5) {

  # Checks
  if (data.table::is.data.table(counts))
    counts <- as.matrix(counts)
  if (!is.matrix(counts))
    counts <- as.matrix(counts)

  if (!is.numeric(response))
    stop("Response must be numeric for regression.")
  if (length(response) != nrow(counts))
    stop("The number of rows in 'counts' must match the length of 'response'.")
  if (CV <= 0 || CV >= 1)
    stop("'CV' must be a value between 0 and 1.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1)
    stop("'alpha' must be a single numeric value in [0, 1].")
  if (!is.numeric(repeats) || repeats < 1)
    stop("'repeats' must be an integer >= 1.")
  repeats <- as.integer(repeats)

  if (is.null(names))
    names <- colnames(counts)
  if (!is.character(names) && !is.factor(names))
    stop("names should be a vector of characters or factors.")
  if (length(names) != ncol(counts))
    stop("length(names) must match ncol(counts).")

  set.seed(seed)

  # Lambda grid
  lambdas <- 10^seq(2, -3, by = -0.1)

  # Output containers
  perf <- vector("list", repeats)
  all_test_sets <- vector("list", repeats)
  all_preds <- vector("list", repeats)
  models <- vector("list", repeats)

  # For each repeat
  for (r in seq_len(repeats)) {

    # Split (expected proportion CV)
    set.seed(seed + r)
    test.set <- sample(c(TRUE, FALSE),
                       size = length(response),
                       replace = TRUE,
                       prob = c(CV, 1 - CV))

    testY <- response[test.set]
    testX <- counts[test.set, , drop = FALSE]
    trainY <- response[!test.set]
    trainX <- counts[!test.set, , drop = FALSE]

    # Inner CV to select lambda
    cvfit <- glmnet::cv.glmnet(
      x = trainX,
      y = trainY,
      alpha = alpha,
      lambda = lambdas,
      standardize = TRUE,
      nfolds = nfolds,
      family = "gaussian"
    )

    lambda_best <- if (use_1se) cvfit$lambda.1se else cvfit$lambda.min

    # Fit final model at chosen lambda
    model <- glmnet::glmnet(
      x = trainX,
      y = trainY,
      alpha = alpha,
      lambda = lambda_best,
      standardize = TRUE,
      family = "gaussian"
    )

    # Predict
    predict_test <- as.vector(predict(model, s = lambda_best, newx = testX))

    test_set_dt <- data.table::data.table(
      repeat_id = r,
      observed = testY,
      predicted = predict_test
    )

    PCC <- stats::cor(test_set_dt$observed, test_set_dt$predicted)
    R2 <- PCC^2

    # Non-zero coefficients ----
    preds <- as.matrix(model$beta)
    preds <- data.table::as.data.table(preds, keep.rownames = "motif")

    # Add display names ----
    add.names <- data.table::data.table(
      motif = colnames(counts),
      name = as.character(names)
    )
    preds <- merge(preds, add.names, by = "motif", all.x = TRUE)
    data.table::setorderv(preds, "s0")
    preds[, repeat_id := r]
    setcolorder(preds, "repeat_id")

    # Return ----
    perf[[r]] <- data.table::data.table(
      repeat_id = r,
      n_test = length(testY),
      n_train = length(trainY),
      PCC = PCC,
      R2 = R2,
      lambda_best = lambda_best
    )
    all_test_sets[[r]] <- test_set_dt
    all_preds[[r]] <- preds
    models[[r]] <- model
  }

  # Aggregate results ----
  perf_all <- data.table::rbindlist(perf)
  test_set_all <- data.table::rbindlist(all_test_sets)
  preds_all <- data.table::rbindlist(all_preds, fill = TRUE)

  # Stability summary ----
  stability <- preds_all[, .(
    selected_n = sum(s0 != 0),
    selected_frac = mean(s0 != 0),
    mean_coef = mean(s0),
    median_coef = stats::median(s0),
    sign_consistency = max(mean(s0 > 0), mean(s0 < 0))
  ), by = .(motif, name)]
  # Best on top
  stability <- stability[order(-selected_frac, -abs(mean_coef))]

  # Return ----
  return(
    list(
      stability = stability,
      performance = perf_all,
      test_set = test_set_all,
      best_predictors = preds_all,
      models = data.table(repeat_id= seq(repeats), models= models),
      alpha = alpha,
      use_1se = use_1se
    )
  )
}
