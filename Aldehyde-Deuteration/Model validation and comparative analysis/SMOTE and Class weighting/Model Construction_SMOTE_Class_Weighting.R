"
# Install
remotes::install_github('https://github.com/barkais/rxn.cond.class.git')
"

library("rxn.cond.class")

if (!requireNamespace("caret", quietly = TRUE)) stop("Please install 'caret' (used for the 75/25 split).")
if (!requireNamespace("knitr", quietly = TRUE)) stop("Please install 'knitr' (used for summary tables).")

n_split_stability <- 10L
stability_seeds <- 1001L:(1000L + n_split_stability)
final_split_seed <- 1011L

pred_accuracy <- function(model, data, resp = "class") {
  p <- predict(model, newdata = data, type = "class")
  mean(as.character(p) == as.character(data[[resp]]), na.rm = TRUE)
}

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

balance_method <- "smote"  # "smote" or "class_weight"
if (!balance_method %in% c("smote", "class_weight")) {
  stop("balance_method must be 'smote' or 'class_weight'")
}

training_csv <- file.path(script_dir, "Training_Data.csv")
external_csv <- file.path(script_dir, "External_Validation_Data.csv")
predict_csv <- file.path(script_dir, "Predicting_New_Substrates_Data.csv")

fit_polr_weighted <- function(formula, data, weights) {
  if (length(weights) != nrow(data)) stop("weights must have length nrow(data)")
  f <- as.formula(formula)
  num.of.vars <- stringi::stri_count(formula, fixed = "+")
  start <- c(rep(0, num.of.vars + 2), 1)
  success <- FALSE
  while (!success) {
    res <- tryCatch(
      list(
        ok = TRUE,
        model = MASS::polr(
          f,
          data = data,
          weights = weights,
          Hess = TRUE,
          start = start,
          control = list(maxit = 100)
        )
      ),
      error = function(e) list(ok = FALSE, err = e)
    )
    if (isTRUE(res$ok)) {
      success <- TRUE
      return(res$model)
    }
    start <- c(0, start)
  }
}

make_class_weights_vec <- function(y) {
  y <- as.factor(y)
  tab <- table(y)
  K <- nlevels(y)
  n <- length(y)
  fac <- as.character(y)
  as.numeric(n / (K * as.numeric(tab[fac])))
}


apply_smote_train <- function(train_data, resp = "class") {
  if (!requireNamespace("recipes", quietly = TRUE) || !requireNamespace("themis", quietly = TRUE)) {
    stop("For balance_method='smote', install.packages(c('recipes', 'themis'))")
  }
  resp <- gsub("`", "", resp, fixed = TRUE)
  if (!resp %in% names(train_data)) stop("response column not found in train_data")
  pred_cols <- setdiff(names(train_data), resp)
  if (!length(pred_cols)) stop("SMOTE requires at least one predictor column besides the response")
  bq <- function(nm) paste0("`", gsub("`", "", nm, fixed = TRUE), "`")
  f <- stats::reformulate(termlabels = vapply(pred_cols, bq, character(1), USE.NAMES = FALSE), response = bq(resp))
  rec0 <- recipes::recipe(f, data = train_data)
  rec1 <- themis::step_smote(rec0, recipes::all_outcomes(), over_ratio = 1)
  rec2 <- recipes::prep(rec1, training = train_data)
  recipes::bake(rec2, new_data = NULL)
}

# Deuteration (SMOTE / class-weighting model construction) -----------------

data <- data.frame(data.table::fread(training_csv), check.names = F)

row.names(data) <- data[,2]
data$class <- as.factor(data$class)
data <- data[,-c(1:2)]

log_msg("=== Setup ===")
cat(sprintf("balance_method: %s\n", balance_method))
cat(sprintf("Samples: n_total=%d (full training pool)\n", nrow(data)))
cat("Class counts (full pool):\n")
print(table(data$class))

# One-off illustration of SMOTE output (same seed / 75%% split as iteration 1) before the stability loop.
if (balance_method == "smote") {
  log_msg(sprintf(
    "=== Example: data after SMOTE (75%% train, seed %d) before stability iterations ===",
    stability_seeds[[1L]]
  ))
  set.seed(stability_seeds[[1L]])
  ex_train_idx <- caret::createDataPartition(data$class, p = 0.75, list = FALSE)
  Train.ex <- data[ex_train_idx, , drop = FALSE]
  cat(sprintf("Rows before SMOTE (illustrative train split): %d\n", nrow(Train.ex)))
  cat("Class counts (before SMOTE):\n")
  print(table(Train.ex$class))
  Train.ex.sm <- tryCatch(
    apply_smote_train(Train.ex),
    error = function(e) {
      warning("SMOTE preview failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(Train.ex.sm)) {
    Train.ex.sm$class <- factor(Train.ex.sm$class, levels = levels(Train.ex$class))
    cat(sprintf("Rows after SMOTE: %d\n", nrow(Train.ex.sm)))
    cat("Class counts (after SMOTE):\n")
    print(table(Train.ex.sm$class))
    cat("First 6 rows after SMOTE (includes synthetic minority rows; row names may be generic):\n")
    print(utils::head(Train.ex.sm, 6L))
  }
} else {
  cat("balance_method is not \"smote\" — no SMOTE preview (class_weight uses original row weights only).\n")
}

# --- Stability across n_split_stability random 75/25 splits (RandomSampling_Final.rmd style) ---
log_msg(sprintf(
  "=== Model search stability: %d stratified splits (seeds %d-%d), top-10 per split ===",
  n_split_stability, min(stability_seeds), max(stability_seeds)
))

pb <- utils::txtProgressBar(
  min = 0,
  max = n_split_stability,
  style = 3,
  width = 50L,
  char = "="
)

split_top_rows <- list()
for (si in seq_along(stability_seeds)) {
  set.seed(stability_seeds[[si]])
  train_idx <- caret::createDataPartition(data$class, p = 0.75, list = FALSE)
  Train.s <- data[train_idx, , drop = FALSE]

  Train.for.s <- Train.s
  if (balance_method == "smote") {
    Train.for.s <- tryCatch(
      apply_smote_train(Train.s),
      error = function(e) {
        warning("SMOTE failed (split ", si, "); using unbalanced train. ", conditionMessage(e))
        Train.s
      }
    )
    Train.for.s$class <- factor(Train.for.s$class, levels = levels(Train.s$class))
  }

  w_sub <- if (balance_method == "class_weight") make_class_weights_vec(Train.s$class) else NULL

  models.s <- sub_model_log(
    data = Train.for.s,
    min = 4,
    max = 4,
    ordinal = TRUE,
    weights = w_sub
  )
  mcf_name <- names(models.s)[grep("McFadden", names(models.s), ignore.case = TRUE)][1L]
  if (is.na(mcf_name) || !nzchar(mcf_name)) {
    stop("sub_model_log output must contain a 'McFadden R2' column.")
  }
  split_top_rows[[si]] <- data.frame(
    split = si,
    seed = stability_seeds[[si]],
    formula = models.s$formula,
    McFadden_R2 = models.s[[mcf_name]],
    stringsAsFactors = FALSE
  )
  utils::setTxtProgressBar(pb, si)
}

close(pb)
cat("\n")

all_top <- do.call(rbind, split_top_rows)
rownames(all_top) <- NULL

form_chr <- as.character(all_top$formula)
uform <- unique(form_chr)
count_per <- vapply(uform, function(f) sum(form_chr == f), integer(1))
avg_mcf <- vapply(uform, function(f) mean(all_top$McFadden_R2[form_chr == f], na.rm = TRUE), numeric(1))
sd_mcf <- vapply(uform, function(f) stats::sd(all_top$McFadden_R2[form_chr == f], na.rm = TRUE), numeric(1))
summary_table <- data.frame(
  formula = uform,
  Count = as.integer(count_per),
  Avg_McFadden = round(avg_mcf, 3),
  SD_McFadden = round(sd_mcf, 3),
  stringsAsFactors = FALSE
)
summary_table <- summary_table[order(-summary_table$Count, -summary_table$Avg_McFadden), ]
rownames(summary_table) <- NULL
summary_table <- head(summary_table, 15L)

log_msg("=== Summary: formulas in top-10 across splits (Count = appearances; tie-break Avg_McFadden) ===")
print(knitr::kable(
  summary_table,
  caption = "Stability summary (cf. RandomSampling_Final.rmd): sort by Count then Avg_McFadden",
  row.names = FALSE
))

test.form <- summary_table$formula[1L]
log_msg(sprintf("Chosen formula (row 1): %s", test.form))

# --- Final 75/25 split for fitting & reporting (same seed convention as Rmd chosen-model chunk) ---
set.seed(final_split_seed)
train_idx <- caret::createDataPartition(data$class, p = 0.75, list = FALSE)
Train.data <- data[train_idx, , drop = FALSE]
Test.data <- data[-train_idx, , drop = FALSE]

cat(sprintf("\nFinal split seed: %d | train n=%d | holdout n=%d\n", final_split_seed, nrow(Train.data), nrow(Test.data)))
cat("Holdout row names:\n")
cat(paste(rownames(Test.data), collapse = ", "), "\n")
cat("Class counts (75%% train, final split):\n")
print(table(Train.data$class))
cat("Class counts (25%% holdout, final split):\n")
print(table(Test.data$class))

Train.for.modeling <- Train.data
if (balance_method == "smote") {
  Train.for.modeling <- tryCatch(
    apply_smote_train(Train.data),
    error = function(e) {
      warning("SMOTE failed on final split; fitting on unbalanced training data. ", conditionMessage(e))
      Train.data
    }
  )
  Train.for.modeling$class <- factor(Train.for.modeling$class, levels = levels(Train.data$class))
  log_msg("SMOTE on final 75% train (all features) before polr fit")
  cat("Class counts on modeling table after SMOTE:\n")
  print(table(Train.for.modeling$class))
  cat(sprintf("Total rows: %d (original train: %d)\n", nrow(Train.for.modeling), nrow(Train.data)))
}

if (balance_method == "class_weight") {
  w_by_class <- tapply(make_class_weights_vec(Train.data$class), Train.data$class, function(x) x[[1L]])
  cat("Case weights applied per class (final split; n/(K*n_c)):\n")
  print(round(w_by_class, 6))
}

if (balance_method == "smote") {
  test <- fit_polr(formula = test.form, data = Train.for.modeling)
} else {
  w <- make_class_weights_vec(Train.data$class)
  test <- fit_polr_weighted(test.form, Train.data, w)
}

log_msg("=== Fitted model (stability-chosen formula on final split) ===")
cat(sprintf("Selected formula: %s\n", test.form))
if (balance_method == "smote") {
  cat(sprintf(
    "Training rows used for fit (after SMOTE): %d (original train rows: %d)\n",
    nrow(Train.for.modeling), nrow(Train.data)
  ))
}
print(summary(test))

# --------- #
# --Train-- #
# --------- #

model.info <- mod.info(test, Train.data, FALSE, FALSE)

log_msg("=== Performance: training set (final split, 75% train) ===")
cat(sprintf("Accuracy: %.2f%%\n", pred_accuracy(test, Train.data) * 100))
cat("Confusion table:\n")
print(model.info$class.table)

suppressWarnings({
  ct_plot(model.info$class.table,
          plot.title = sprintf("Training set (final split, seed %d)", final_split_seed),
          conformation = '1. 1st Place')$plot

  prob.heatmap(test, Train.data,
               plot.title = sprintf("Training set (final split, seed %d)", final_split_seed),
               conformation = '1. 1st Place')
})

# ---------- #
# -- Test -- #
# ---------- #

model.info <- mod.info(test, Test.data, FALSE, FALSE)

log_msg("=== Performance: holdout test (final split, 25%) ===")
cat(sprintf("Accuracy: %.2f%%\n", pred_accuracy(test, Test.data) * 100))
cat("Confusion table:\n")
print(model.info$class.table)

suppressWarnings({
  ct_plot(model.info$class.table,
          plot.title = sprintf("Test set (final split, seed %d)", final_split_seed),
          conformation = '1. 1st Place')$plot

  prob.heatmap(test, Test.data,
               plot.title = sprintf("Test set (final split, seed %d)", final_split_seed),
               conformation = '1. 1st Place')
})

# ------------------------- #
# -- External Validation -- #
# ------------------------- #

External <- data.frame(data.table::fread(external_csv), check.names = F)
RN <- External$V1
External <- External[,-1]
External$class <- as.factor(External$class)
row.names(External) <- RN
colnames(External) <- colnames(Train.data[, c(2:dim(Train.data)[2])])

model.info <- mod.info(test, External, FALSE, FALSE)

log_msg("=== Performance: external validation ===")
cat(sprintf("Accuracy: %.2f%%\n", pred_accuracy(test, External) * 100))
cat("Confusion table:\n")
print(model.info$class.table)

suppressWarnings({
  ct_plot(model.info$class.table,
          plot.title = 'External Validation',
          conformation = '1. 1st Place')$plot

  prob.heatmap(test, External,
               plot.title = 'External Validation',
               conformation = '1. 1st Place')
})

# ----------------- #
# -- Predictions -- #
# ----------------- #

Prediction.set <- data.frame(data.table::fread(predict_csv), check.names = F)
RN <- Prediction.set$V1
Prediction.set <- Prediction.set[,-1]
row.names(Prediction.set) <- RN

log_msg("=== Predictions (new substrates) ===")
pred_tbl <- cbind(
  predict(test, Prediction.set, "probs") * 100,
  predicted_class = predict(test, Prediction.set, "class")
)
print(knitr::kable(pred_tbl, digits = 3))
