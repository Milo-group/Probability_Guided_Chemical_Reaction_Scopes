"
# Install
remotes::install_github('https://github.com/barkais/rxn.cond.class.git')
"

library("rxn.cond.class")

if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("Please install 'knitr' (used for summary tables).")
}

script_dir <- dirname(normalizePath(if (is.null(sys.frame(1)$ofile)) "." else sys.frame(1)$ofile))
source(file.path(script_dir, "kennard_stone_sampling.R"))

training_csv <- file.path(script_dir, "Training_Data.csv")
external_csv <- file.path(script_dir, "External_Validation_Data.csv")
predict_csv <- file.path(script_dir, "Predicting_New_Substrates_Data.csv")

# Per-class Kennard–Stone picks (up to n_per_class; smaller classes keep all rows).
n_per_class <- NULL  # NULL = smallest class size (cf. default simi.sampler coverage)
top_k_models <- 10L

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

pred_accuracy <- function(model, data, resp = "class") {
  p <- predict(model, newdata = data, type = "class")
  mean(as.character(p) == as.character(data[[resp]]), na.rm = TRUE)
}

mcfadden_col_name <- function(models_df) {
  nm <- names(models_df)[grep("McFadden", names(models_df), ignore.case = TRUE)][1L]
  if (is.na(nm) || !nzchar(nm)) {
    stop("sub_model_log output must contain a 'McFadden R2' column.", call. = FALSE)
  }
  nm
}

calc_mcfadden_r2 <- function(model, eval_data, response_col = "class") {
  y <- factor(eval_data[[response_col]])
  probs <- predict(model, newdata = eval_data, type = "probs")

  if (!is.matrix(probs)) {
    probs <- as.matrix(probs)
    if (ncol(probs) == 1L) {
      colnames(probs) <- levels(y)
    }
  }

  probs <- probs[, levels(y), drop = FALSE]
  row_id <- seq_len(nrow(eval_data))
  col_id <- match(as.character(y), colnames(probs))
  p_full <- pmax(probs[cbind(row_id, col_id)], .Machine$double.eps)
  ll_full <- sum(log(p_full))

  null_probs <- as.numeric(prop.table(table(y)))
  names(null_probs) <- levels(y)
  p_null <- pmax(null_probs[as.character(y)], .Machine$double.eps)
  ll_null <- sum(log(p_null))

  1 - (ll_full / ll_null)
}

fit_top_models_with_metrics <- function(models_df, train_data, test_data, top_k = 10L) {
  top_k <- min(as.integer(top_k), nrow(models_df))
  if (top_k <= 0L) {
    return(models_df[0, , drop = FALSE])
  }

  train_acc <- numeric(top_k)
  test_acc <- numeric(top_k)
  q2_mcfadden <- numeric(top_k)

  for (i in seq_len(top_k)) {
    form_i <- models_df[i, 1][[1]]
    model_i <- tryCatch(
      fit_polr(formula = form_i, data = train_data),
      error = function(e) NULL
    )
    if (is.null(model_i)) {
      train_acc[i] <- NA_real_
      test_acc[i] <- NA_real_
      q2_mcfadden[i] <- NA_real_
    } else {
      train_acc[i] <- pred_accuracy(model_i, train_data)
      test_acc[i] <- pred_accuracy(model_i, test_data)
      q2_mcfadden[i] <- calc_mcfadden_r2(model_i, test_data, "class")
    }
  }

  out <- models_df[seq_len(top_k), , drop = FALSE]
  out$Train_Accuracy_Pct <- round(train_acc * 100, 2)
  out$Test_Accuracy_Pct <- round(test_acc * 100, 2)
  out$`Q^2_McFadden_Test` <- round(q2_mcfadden, 3)
  out
}

print_model_metrics <- function(model, train_data, test_data, external_data, label_train, label_test) {
  mi_tr <- mod.info(model, train_data, FALSE, FALSE)
  log_msg(sprintf("=== Performance: %s ===", label_train))
  cat(sprintf("Accuracy: %.2f%%\n", pred_accuracy(model, train_data) * 100))
  cat("Confusion table:\n")
  print(mi_tr$class.table)

  mi_te <- mod.info(model, test_data, FALSE, FALSE)
  log_msg(sprintf("=== Performance: %s ===", label_test))
  cat(sprintf("Accuracy: %.2f%%\n", pred_accuracy(model, test_data) * 100))
  cat(sprintf("Q^2 (McFadden on test): %.4f\n", calc_mcfadden_r2(model, test_data, "class")))
  cat("Confusion table:\n")
  print(mi_te$class.table)

  if (!is.null(external_data) && nrow(external_data) > 0) {
    mi_ex <- mod.info(model, external_data, FALSE, FALSE)
    log_msg("=== Performance: external validation ===")
    cat(sprintf("Accuracy: %.2f%%\n", pred_accuracy(model, external_data) * 100))
    cat("Confusion table:\n")
    print(mi_ex$class.table)
  }

  invisible(list(train = mi_tr, test = mi_te))
}

# Deuteration (Kennard–Stone train / holdout split) ---------------------------------

data <- prepare_deuteration_training_data(training_csv)

log_msg("=== Setup (Kennard–Stone sampling) ===")
cat(sprintf("Samples: n_total=%d\n", nrow(data)))
cat("Class counts (full pool):\n")
print(table(data$class))

train_idx <- build_kennard_stone_train_indices(
  data = data,
  class_col = "class",
  n_per_class = n_per_class
)

Train.data <- data[train_idx, , drop = FALSE]
Test.data <- data[-train_idx, , drop = FALSE]

cat(sprintf(
  "Kennard–Stone train n=%d | holdout n=%d (n_per_class=%s)\n",
  nrow(Train.data),
  nrow(Test.data),
  if (is.null(n_per_class)) as.character(min(table(data$class))) else as.character(n_per_class)
))
cat("Class counts (train):\n")
print(table(Train.data$class))
cat("Class counts (holdout):\n")
print(table(Test.data$class))
cat("Holdout row names:\n")
cat(paste(rownames(Test.data), collapse = ", "), "\n")

log_msg(sprintf("=== Model search on Kennard–Stone train set (top %d by McFadden R²) ===", top_k_models))
models.ordinal <- sub_model_log(
  data = Train.data,
  min = 4,
  max = 4,
  ordinal = TRUE
)

if (is.null(models.ordinal) || nrow(models.ordinal) == 0) {
  stop("sub_model_log returned no candidate models on the Kennard–Stone training set.")
}

top_models <- fit_top_models_with_metrics(
  models_df = models.ordinal,
  train_data = Train.data,
  test_data = Test.data,
  top_k = top_k_models
)

log_msg("=== Top candidate models (ranked by sub_model_log / McFadden R²) ===")
print(knitr::kable(top_models, row.names = FALSE))

test.form <- as.character(top_models[1, 1][[1]])
mcf_col <- mcfadden_col_name(top_models)
log_msg(sprintf(
  "Chosen formula (rank 1): %s | McFadden R² = %s | train acc = %s%% | test acc = %s%% | Q^2 = %s",
  test.form,
  format(top_models[[mcf_col]][1], digits = 4),
  format(top_models$Train_Accuracy_Pct[1], digits = 4),
  format(top_models$Test_Accuracy_Pct[1], digits = 4),
  format(top_models$`Q^2_McFadden_Test`[1], digits = 4)
))

test <- fit_polr(formula = test.form, data = Train.data)

log_msg("=== Fitted model (chosen formula on Kennard–Stone train set) ===")
cat(sprintf("Selected formula: %s\n", test.form))
print(summary(test))

External <- data.frame(data.table::fread(external_csv), check.names = FALSE)
RN <- External$V1
External <- External[, -1]
External$class <- as.factor(External$class)
row.names(External) <- RN
colnames(External) <- colnames(Train.data[, c(2:dim(Train.data)[2]), drop = FALSE])

print_model_metrics(
  model = test,
  train_data = Train.data,
  test_data = Test.data,
  external_data = External,
  label_train = "training set (Kennard–Stone)",
  label_test = "holdout test (Kennard–Stone)"
)

# --------- #
# --Train-- #
# --------- #

model.info <- mod.info(test, Train.data, FALSE, FALSE)

ct_plot(
  model.info$class.table,
  plot.title = "Training Set (Kennard–Stone)",
  conformation = "1. 1st Place"
)$plot

prob.heatmap(
  test,
  Train.data,
  plot.title = "Training Set (Kennard–Stone)",
  conformation = "1. 1st Place"
)

# ---------- #
# -- Test -- #
# ---------- #

model.info <- mod.info(test, Test.data, FALSE, FALSE)

ct_plot(
  model.info$class.table,
  plot.title = "Test Set (Kennard–Stone holdout)",
  conformation = "1. 1st Place"
)$plot

prob.heatmap(
  test,
  Test.data,
  plot.title = "Test Set (Kennard–Stone holdout)",
  conformation = "1. 1st Place"
)

# ------------------------- #
# -- External Validation -- #
# ------------------------- #

model.info <- mod.info(test, External, FALSE, FALSE)

ct_plot(
  model.info$class.table,
  plot.title = "External Validation",
  conformation = "1. 1st Place"
)$plot

prob.heatmap(
  test,
  External,
  plot.title = "External Validation",
  conformation = "1. 1st Place"
)

# ----------------- #
# -- Predictions -- #
# ----------------- #

Prediction.set <- data.frame(data.table::fread(predict_csv), check.names = FALSE)
RN <- Prediction.set$V1
Prediction.set <- Prediction.set[, -1]
row.names(Prediction.set) <- RN

log_msg("=== Predictions (new substrates) ===")
pred_tbl <- cbind(
  predict(test, Prediction.set, "probs") * 100,
  predicted_class = predict(test, Prediction.set, "class")
)
print(knitr::kable(pred_tbl, digits = 3))

log_msg("Finished Model Construction_Kennard_Stone.R.")
