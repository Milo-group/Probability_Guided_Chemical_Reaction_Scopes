"
# Install
remotes::install_github('https://github.com/barkais/rxn.cond.class.git')
"

library("rxn.cond.class")

if (!requireNamespace("nnet", quietly = TRUE)) {
  stop("Please install 'nnet' (used for multinomial models).")
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
    probs <- cbind(probs)
    colnames(probs) <- levels(y)
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

fit_top_models_with_metrics <- function(models_df, train_data, test_data, top_k = 15L) {
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
      nnet::multinom(form_i, data = train_data, maxit = 2000, trace = FALSE),
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

print_model_table <- function(tbl) {
  if (requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(tbl, row.names = FALSE))
  } else {
    print(tbl, row.names = FALSE)
  }
}

# Folder with this script and Training_morgan_filtered.csv (not necessarily getwd())
.morgan_features_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) {
    return(dirname(normalizePath(sub("^--file=", "", f[1L]), winslash = "/", mustWork = FALSE)))
  }
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of) && nzchar(of)) {
      return(dirname(normalizePath(of, winslash = "/", mustWork = FALSE)))
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) {
      return(dirname(normalizePath(p, winslash = "/", mustWork = FALSE)))
    }
  }
  getwd()
}

# Hydrogen Isotope Exchange — Morgan FP (redundant bits removed)
# Regenerate CSVs: python morgan_fp_filter.py

.morgan_dir <- .morgan_features_dir()
source(file.path(.morgan_dir, "morgan_feature_config.R"), local = TRUE)
.morgan_load <- load_morgan_feature_table(.morgan_dir)
.morgan_spec <- .morgan_load$spec
training_csv <- .morgan_load$training_csv
prediction_csv <- .morgan_load$prediction_csv
feat_cols <- .morgan_load$feature_cols
data50 <- .morgan_load$data

## Sampler data (same compounds as training; similarity in feature space)
sampl.dat <- data50[, c(feat_cols, "class", "flag"), drop = FALSE]
sampl.dat$class <- factor(sampl.dat$class)

# Sample from group 1 using simi.sampler, returning samples from the smallest group
one <- simi.sampler(sampl.dat, 1)
two <- simi.sampler(sampl.dat, 2)
three <- simi.sampler(sampl.dat, 3)
four <- simi.sampler(sampl.dat, 4)

# Sample from group 2 molecules that are similar to group 1
two_three <- simi.sampler(sampl.dat, 2, 3)

# Sample from group 1 molecules that are similar to group 3
one_three <- simi.sampler(sampl.dat, 1, 3)

# Sample from group 4 molecules that are similar to group 3
four_three <- simi.sampler(sampl.dat, 4, 3)

print(one)
print(two)
print(three)
print(four)
print(two_three)
print(one_three)
print(four_three)

## Back to model search

similarties <- c(
  three,
  union(two, two_three),
  union(one, one_three),
  union(four, four_three)
)

Train.data <- data50[similarties, , drop = FALSE]
Test.data <- data50[-similarties, , drop = FALSE]

# Non-ordinal multinomial logistic regression models (subset size from feature method)
models <- sub_model_log(
  data = Train.data,
  min = .morgan_spec$model_min,
  max = .morgan_spec$model_max,
  ordinal = FALSE
)

top_k_models <- 15L
top_models <- fit_top_models_with_metrics(
  models_df = models,
  train_data = Train.data,
  test_data = Test.data,
  top_k = top_k_models
)

test.form <- top_models[3, 1][[1]]
mcf_col <- mcfadden_col_name(top_models)

test <- nnet::multinom(
  test.form,
  data = Train.data,
  maxit = 2000,
  trace = FALSE
)

cat("\n", strrep("=", 64), "\n", "FITTED MODEL (", .morgan_spec$label, " + multinomial)\n", strrep("=", 64), "\n", sep = "")
cat("\nFormula (chosen row 3 of sub_model_log ranking):\n")
print(stats::as.formula(if (inherits(test.form, "formula")) test.form else paste(test.form)))
cat(sprintf(
  "\nChosen model metrics: McFadden R2 (train search) = %s | train acc = %s%% | test acc = %s%% | Q^2 (McFadden test) = %s\n",
  format(top_models[[mcf_col]][3], digits = 4),
  format(top_models$Train_Accuracy_Pct[3], digits = 4),
  format(top_models$Test_Accuracy_Pct[3], digits = 4),
  format(top_models$`Q^2_McFadden_Test`[3], digits = 4)
))
cat("\nMcFadden R2 ranking — top ", top_k_models, " candidate formulas (train/test accuracy and Q^2 on holdout):\n", sep = "")
print_model_table(top_models)
cat("\n--- summary(test) ---\n")
print(summary(test))
cat("\nResidual deviance:", stats::deviance(test), "  AIC:", stats::AIC(test), "\n")

# --------- #
# --Train-- #
# --------- #

model.info.train <- mod.info(test, Train.data, FALSE, FALSE)
cat("\n", strrep("=", 64), "\n", "TRAINING SET PERFORMANCE (n = ", nrow(Train.data), ")\n", strrep("=", 64), "\n", sep = "")
cat(sprintf("Accuracy: %.4f (%.2f%%)\n", model.info.train$accuracy, model.info.train$accuracy * 100))
if (!is.null(model.info.train$McFadden)) {
  cat(sprintf("McFadden R2 (vs. training null): %.4f\n", model.info.train$McFadden))
}
cat(sprintf("McFadden R2 (probability-based): %.4f\n", calc_mcfadden_r2(test, Train.data, "class")))
cat("\nConfusion matrix (rows = observed class, cols = predicted):\n")
print(model.info.train$class.table)

ct_plot(
  model.info.train$class.table,
  plot.title = "Training Set",
  conformation = .morgan_spec$label
)$plot

prob.heatmap(
  test,
  Train.data,
  plot.title = "Training Set",
  conformation = .morgan_spec$label
)

# ---------- #
# -- Test -- #
# ---------- #

model.info.test <- mod.info(test, Test.data, FALSE, FALSE)
cat("\n", strrep("=", 64), "\n", "HOLDOUT TEST SET PERFORMANCE (n = ", nrow(Test.data), ")\n", strrep("=", 64), "\n", sep = "")
cat(sprintf("Accuracy: %.4f (%.2f%%)\n", model.info.test$accuracy, model.info.test$accuracy * 100))
if (!is.null(model.info.test$McFadden)) {
  cat(sprintf("McFadden R2: %.4f\n", model.info.test$McFadden))
}
cat(sprintf("Q^2 (McFadden on test): %.4f\n", calc_mcfadden_r2(test, Test.data, "class")))
cat("\nConfusion matrix (rows = observed class, cols = predicted):\n")
print(model.info.test$class.table)

ct_plot(
  model.info.test$class.table,
  plot.title = "Test Set",
  conformation = .morgan_spec$label
)$plot

prob.heatmap(
  test,
  Test.data,
  plot.title = "Test Set",
  conformation = .morgan_spec$label
)

# ----------------- #
# -- Predictions -- #
# ----------------- #

if (!file.exists(prediction_csv)) {
  message("Missing ", prediction_csv, " — skipping prediction table.")
  cat("\n", strrep("=", 64), "\nPREDICTIONS: skipped (file not found)\n", strrep("=", 64), "\n", sep = "")
} else {
  pred_raw <- data.frame(data.table::fread(prediction_csv), check.names = FALSE)
  miss_feat <- setdiff(feat_cols, names(pred_raw))
  if (length(miss_feat)) {
    stop(prediction_csv, " must contain the same Morgan FP columns as ", training_csv)
  }
  row.names(pred_raw) <- as.character(pred_raw$Name)
  pred_raw$Name <- NULL
  Prediction.set <- pred_raw[, feat_cols, drop = FALSE]
  for (j in feat_cols) Prediction.set[[j]] <- as.numeric(Prediction.set[[j]])

  cls <- predict(test, Prediction.set, type = "class")
  cat("\n", strrep("=", 64), "\n", "PREDICTED CLASSES (n = ", nrow(Prediction.set), ")\n", strrep("=", 64), "\n", sep = "")
  print(data.frame(
    Name = row.names(Prediction.set),
    predicted_class = as.character(cls),
    stringsAsFactors = FALSE
  ), row.names = FALSE)
}
