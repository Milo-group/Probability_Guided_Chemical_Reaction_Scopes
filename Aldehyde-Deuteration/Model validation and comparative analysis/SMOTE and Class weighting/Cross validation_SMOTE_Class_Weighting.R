library('rxn.cond.class')

if (!requireNamespace("caret", quietly = TRUE)) stop("Please install 'caret' (used for the 75/25 split and CV folds).")

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

pred_accuracy <- function(model, data, resp = "class") {
  p <- predict(model, newdata = data, type = "class")
  mean(as.character(p) == as.character(data[[resp]]), na.rm = TRUE)
}

balance_method <- "smote"  # "smote" or "class_weight"
if (!balance_method %in% c("smote", "class_weight")) {
  stop("balance_method must be 'smote' or 'class_weight'")
}


script_dir <- dirname(normalizePath(if (is.null(sys.frame(1)$ofile)) "." else sys.frame(1)$ofile))
training_csv <- file.path(script_dir, "Training_Data.csv")

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

# SMOTE on in-fold training rows using every predictor column (immediately after the fold
# train split, before polr / multinom fit). Held-out fold rows are never SMOTE'd.
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

# Cross validation code ---------------------------------------------------

cv_classification <- function(data,
                              test.form,
                              ordinal = FALSE,
                              k = 5,
                              n.iter = 1,
                              seed = NULL,
                              verbose = FALSE,
                              balance_method = "class_weight",
                              ...) {
  if (!requireNamespace("caret", quietly = TRUE)) stop("Please install 'caret'")
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Please install 'nnet'")
  
  test.form1 <- as.formula(test.form)
  
  resp_var <- deparse(test.form1[[2]])
  resp_var_clean <- gsub("`", "", resp_var)
  
  if (!resp_var_clean %in% names(data)) stop("response variable not found in data")
  if (!is.factor(data[[resp_var_clean]])) {
    warning("Response variable coerced to factor for stratified sampling.")
    data[[resp_var_clean]] <- factor(data[[resp_var_clean]])
  }
  
  n <- nrow(data)
  classes <- levels(data[[resp_var_clean]])
  
  loo <- FALSE
  if (is.character(k) && tolower(k) == "loo") {
    k <- n
    loo <- TRUE
  }
  
  if (k == n && n.iter > 1) {
    if (verbose) message("LOO detected. Repeats ignored; performing LOO once.")
    n.iter <- 1
  }
  
  results_table <- data.frame()
  
  total_folds <- n.iter * k
  pb <- utils::txtProgressBar(min = 0, max = total_folds, style = 3)
  progress_counter <- 0
  
  for (iter in seq_len(n.iter)) {
    if (!is.null(seed)) set.seed(seed + iter)
    
    if (k == n) {
      folds <- as.list(seq_len(n))
    } else {
      folds <- caret::createFolds(data[[resp_var_clean]], k = k, list = TRUE, returnTrain = FALSE)
    }
    
    for (f in seq_along(folds)) {
      test_idx <- folds[[f]]
      Train.data <- data[-test_idx, , drop = FALSE]
      Test.data  <- data[test_idx, , drop = FALSE]
      
      Train.data[[resp_var_clean]] <- factor(Train.data[[resp_var_clean]], levels = classes)
      
      if (length(unique(Train.data[[resp_var_clean]])) < 2) next
      
      # SMOTE in-fold training rows (all features) before fitting; test fold unchanged.
      train_for_fit <- Train.data
      if (ordinal && balance_method == "smote") {
        train_for_fit <- tryCatch(
          apply_smote_train(Train.data, resp = resp_var_clean),
          error = function(e) {
            if (verbose) warning(conditionMessage(e))
            Train.data
          }
        )
        train_for_fit[[resp_var_clean]] <- factor(train_for_fit[[resp_var_clean]], levels = classes)
      }
      
      if (!ordinal) {
        model <- nnet::multinom(test.form1, data = train_for_fit, maxit = 2000, trace = FALSE)
      } else {
        if (balance_method == "class_weight") {
          w <- make_class_weights_vec(train_for_fit[[resp_var_clean]])
          model <- fit_polr_weighted(test.form, train_for_fit, w)
        } else {
          model <- fit_polr(formula = test.form, data = train_for_fit)
        }
      }
      
      pred_class <- predict(model, newdata = Test.data, type = "class")
      acc <- mean(as.character(pred_class) == as.character(Test.data[[resp_var_clean]]))
      
      class_counts <- table(factor(Test.data[[resp_var_clean]], levels = classes))
      
      row <- data.frame(
        iteration = iter,
        fold = f,
        left_out_samples = paste(test_idx, collapse = ","),
        accuracy = acc
      )
      
      if (loo) row$predicted_class <- as.character(pred_class)
      
      row <- cbind(row, as.list(as.numeric(class_counts)))
      names(row)[(ncol(row)-length(classes)+1):ncol(row)] <- classes
      
      results_table <- rbind(results_table, row)
      
      progress_counter <- progress_counter + 1
      utils::setTxtProgressBar(pb, progress_counter)
      
      if (verbose) {
        message(sprintf("  Iter %d Fold %d accuracy = %.4f", iter, f, acc))
      }
    }
  }
  
  close(pb)
  
  overall_mean <- mean(results_table$accuracy, na.rm = TRUE)
  
  list(
    results_table = results_table,
    overall_mean_accuracy = overall_mean,
    k = if (k == n) "LOO" else k,
    n.iter = n.iter
  )
}

# Deuteration (SMOTE / class-weighting validation) ------------------------

data <- data.frame(data.table::fread(training_csv), check.names = F)

row.names(data) <- data[,2]
data$class <- as.factor(data$class)
data <- data[,-c(1:2)]

# CHANGED: random 75% / 25% stratified split replaces `simi.sampler` train/test construction.
set.seed(10)
train_idx <- caret::createDataPartition(data$class, p = 0.75, list = FALSE)
Train.data <- data[train_idx, , drop = FALSE]
Test.data <- data[-train_idx, , drop = FALSE]

test.form <- "`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`"

log_msg("=== Setup ===")
cat(sprintf("balance_method: %s\n", balance_method))
cat(sprintf("Formula: %s\n", test.form))
cat(sprintf(
  "Samples: n_total=%d | train (75%%)=%d | holdout test (25%%)=%d\n",
  nrow(data), nrow(Train.data), nrow(Test.data)
))
cat("Class counts (full data):\n")
print(table(data$class))
cat("Class counts (75%% train):\n")
print(table(Train.data$class))
cat("Class counts (25%% holdout):\n")
print(table(Test.data$class))
if (balance_method == "class_weight") {
  w_setup <- make_class_weights_vec(Train.data$class)
  w_by_class <- tapply(w_setup, Train.data$class, FUN = function(x) x[[1L]])
  cat("Case weights per class for full 75%% train pool (n/(K*n_c)); each CV fold recomputes from in-fold counts:\n")
  print(round(w_by_class, 6))
}

CV7 <- cv_classification(
  Train.data,
  test.form,
  ordinal = TRUE,
  k = 7,
  n.iter = 200,
  seed = 10,
  balance_method = balance_method
)

LOO <- cv_classification(
  Train.data,
  test.form,
  ordinal = TRUE,
  k = "loo",
  n.iter = 200,
  seed = 10,
  balance_method = balance_method
)

log_msg("=== Cross-validation performance (on 75% training pool) ===")
acc7 <- CV7$results_table$accuracy
accl <- LOO$results_table$accuracy
cat(sprintf(
  "7-fold CV: mean accuracy = %.2f%% | SD = %.2f%% | min = %.2f%% | max = %.2f%%\n",
  mean(acc7, na.rm = TRUE) * 100,
  stats::sd(acc7, na.rm = TRUE) * 100,
  min(acc7, na.rm = TRUE) * 100,
  max(acc7, na.rm = TRUE) * 100
))
cat(sprintf(
  "LOO: mean accuracy = %.2f%% | SD = %.2f%% | min = %.2f%% | max = %.2f%%\n",
  mean(accl, na.rm = TRUE) * 100,
  stats::sd(accl, na.rm = TRUE) * 100,
  min(accl, na.rm = TRUE) * 100,
  max(accl, na.rm = TRUE) * 100
))

# -------------------------------
# Model report (PDF) generation
# -------------------------------
report_dir <- file.path(dirname(dirname(script_dir)), "Model Reports")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)
report_file <- file.path(
  report_dir,
  sprintf("CrossValidation.ModelReport.SMOTE_ClassWeighting.%s.pdf", balance_method)
)

log_msg(sprintf("Writing cross-validation model report PDF: %s", report_file))

pdf(report_file, width = 11, height = 8.5)

acc_cv_pct <- CV7$results_table$accuracy * 100
acc_loo_pct <- LOO$results_table$accuracy * 100

par(mar = c(5, 4, 4, 2) + 0.1)
boxplot(
  list("CV (k=7)" = acc_cv_pct, "LOO" = acc_loo_pct),
  main = sprintf("Cross-validation accuracy (%s)", balance_method),
  ylab = "Accuracy (%)",
  col = c("#4C78A8", "#F58518")
)

hist(acc_cv_pct, breaks = 20, col = "#4C78A8AA", border = "#4C78A8",
     main = "CV (k=7) accuracy histogram", xlab = "Accuracy (%)")
hist(acc_loo_pct, breaks = 20, col = "#F58518AA", border = "#F58518",
     main = "LOO accuracy histogram", xlab = "Accuracy (%)")


train_for_final <- Train.data
if (balance_method == "smote") {
  train_for_final <- tryCatch(
    apply_smote_train(Train.data),
    error = function(e) {
      warning("SMOTE failed on full training subset; fitting on unbalanced data. ", conditionMessage(e))
      Train.data
    }
  )
  train_for_final$class <- factor(train_for_final$class, levels = levels(Train.data$class))
  log_msg("Class counts after SMOTE (75% training subset, all features, final model fit)")
  print(table(train_for_final$class))
  cat(sprintf("Total rows after SMOTE: %d (before SMOTE: %d)\n", nrow(train_for_final), nrow(Train.data)))
  final_model <- fit_polr(formula = test.form, data = train_for_final)
} else {
  w <- make_class_weights_vec(Train.data$class)
  final_model <- fit_polr_weighted(test.form, Train.data, w)
}

log_msg("=== Final fitted model (same formula & balance method as CV) ===")
print(summary(final_model))
cat(sprintf(
  "Accuracy on 75%% train (original rows, not SMOTE duplicates): %.2f%%\n",
  pred_accuracy(final_model, Train.data) * 100
))
if (nrow(Test.data) > 0) {
  cat(sprintf("Accuracy on 25%% holdout: %.2f%%\n", pred_accuracy(final_model, Test.data) * 100))
}
mi_tr <- mod.info(final_model, Train.data, FALSE, FALSE)
cat("Confusion table — train (75%%):\n")
print(mi_tr$class.table)
if (nrow(Test.data) > 0) {
  mi_te <- mod.info(final_model, Test.data, FALSE, FALSE)
  cat("Confusion table — holdout test (25%%):\n")
  print(mi_te$class.table)
}

grid_side_by_side_confusion <- function(confusion_plot, confusion_table, left_title, right_title) {
  if (!requireNamespace("grid", quietly = TRUE)) stop("`grid` package is required.")
  
  confusion_grob <- confusion_plot
  if (inherits(confusion_plot, "ggplot")) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("`ggplot2` is required to render ct_plot() output.")
    }
    confusion_grob <- ggplot2::ggplotGrob(confusion_plot)
  }
  
  tab_mat <- as.matrix(confusion_table)
  tab_lines <- capture.output(print(tab_mat))
  tab_label <- paste(tab_lines, collapse = "\n")
  
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(confusion_grob)
  grid::grid.text(left_title, x = 0.5, y = 0.98, just = "top",
                  gp = grid::gpar(fontface = "bold", cex = 1.0))
  grid::popViewport()
  
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.text(right_title, x = 0.5, y = 0.98, just = "top",
                  gp = grid::gpar(fontface = "bold", cex = 1.0))
  grid::grid.text(tab_label, x = 0.02, y = 0.96, just = "left",
                  gp = grid::gpar(fontfamily = "mono", cex = 0.65),
                  default.units = "npc")
  grid::popViewport()
  
  grid::popViewport()
}

cm_train <- ct_plot(
  mi_tr$class.table,
  plot.title = "Train Set (75% random split; for model report)",
  conformation = "1. 1st Place"
)
grid_side_by_side_confusion(
  confusion_plot = cm_train$plot,
  confusion_table = mi_tr$class.table,
  left_title = "Train confusion matrix",
  right_title = "Train confusion table"
)

p_train_hm <- prob.heatmap(
  final_model,
  Train.data,
  plot.title = "Train Set (75% random split; for model report)",
  conformation = "1. 1st Place"
)
print(p_train_hm)

if (nrow(Test.data) > 0) {
  cm_test <- ct_plot(
    mi_te$class.table,
    plot.title = "Test Set (25% holdout; for model report)",
    conformation = "1. 1st Place"
  )

  grid_side_by_side_confusion(
    confusion_plot = cm_test$plot,
    confusion_table = mi_te$class.table,
    left_title = "Test confusion matrix",
    right_title = "Test confusion table"
  )

  p_test_hm <- prob.heatmap(
    final_model,
    Test.data,
    plot.title = "Test Set (25% holdout; for model report)",
    conformation = "1. 1st Place"
  )
  print(p_test_hm)
}

dev.off()
log_msg("Cross-validation model report PDF written successfully.")

log_msg("Finished Cross validation_SMOTE_Class_Weighting.R.")
