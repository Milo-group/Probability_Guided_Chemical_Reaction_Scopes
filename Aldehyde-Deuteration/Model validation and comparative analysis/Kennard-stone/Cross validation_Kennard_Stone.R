library("rxn.cond.class")

if (!requireNamespace("caret", quietly = TRUE)) {
  stop("Please install 'caret' (used for CV folds).")
}

script_dir <- dirname(normalizePath(if (is.null(sys.frame(1)$ofile)) "." else sys.frame(1)$ofile))
source(file.path(script_dir, "kennard_stone_sampling.R"))

training_csv <- file.path(script_dir, "Training_Data.csv")
external_csv <- file.path(script_dir, "External_Validation_Data.csv")

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

pred_accuracy <- function(model, data, resp = "class") {
  p <- predict(model, newdata = data, type = "class")
  mean(as.character(p) == as.character(data[[resp]]), na.rm = TRUE)
}

load_external_validation_aligned <- function(external_csv, ref_data) {
  if (is.null(external_csv) || !nzchar(external_csv) || !file.exists(external_csv)) {
    return(NULL)
  }
  External <- data.frame(data.table::fread(external_csv, check.names = FALSE))
  RN <- External[[1]]
  External <- External[, -1, drop = FALSE]
  row.names(External) <- RN
  n_ref <- ncol(ref_data)
  if (ncol(External) == n_ref - 1L) {
    colnames(External) <- colnames(ref_data)[seq_len(n_ref)[-1L]]
  } else if (ncol(External) == n_ref) {
    colnames(External) <- colnames(ref_data)
  } else {
    warning(sprintf(
      "External_Validation_Data.csv has %d columns after row-id drop; expected %d (no flag) or %d (full). Skipping external validation.",
      ncol(External), n_ref - 1L, n_ref
    ))
    return(NULL)
  }
  External$class <- factor(as.character(External$class), levels = levels(ref_data$class))
  External
}

build_prediction_table <- function(model, data) {
  probs <- predict(model, data, "probs") * 100
  pred_class <- predict(model, data, "class")
  cbind(
    actual_class = as.character(data$class),
    as.data.frame(probs),
    predicted_class = as.character(pred_class)
  )
}

prob_heatmap_scale_params <- function(n_rows) {
  n_rows <- max(as.integer(n_rows), 1L)
  text_size <- max(6, 10 - (n_rows - 10) / 5)
  plot_ratio <- if (n_rows > 10) 0.3 + (10 / n_rows) else 0.5
  list(
    text_size = text_size,
    plot_ratio = plot_ratio,
    cell_text_size = text_size / 2.5,
    exp_text_size = text_size / 2
  )
}

scale_prob_heatmap_plot <- function(p, n_rows) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`ggplot2` is required to scale probability heatmaps.")
  }
  scale <- prob_heatmap_scale_params(n_rows)
  y_col <- p$theme$axis.text.y$colour

  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomText")) {
      layer_size <- p$layers[[i]]$aes_params$size
      p$layers[[i]]$aes_params$size <- if (is.null(layer_size)) {
        scale$cell_text_size
      } else {
        scale$exp_text_size
      }
    }
  }

  title_size <- max(9, min(13, scale$text_size + 2))
  p +
    ggplot2::coord_fixed(ratio = scale$plot_ratio) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = scale$text_size, face = "bold"),
      axis.text.y = ggplot2::element_text(
        size = scale$text_size,
        face = "bold",
        colour = y_col
      ),
      plot.title = ggplot2::element_text(size = title_size),
      plot.subtitle = ggplot2::element_text(size = title_size * 0.85),
      plot.caption = ggplot2::element_text(size = max(7, scale$text_size * 0.75)),
      legend.text = ggplot2::element_text(size = max(5, scale$text_size * 0.65)),
      legend.key.height = ggplot2::unit(max(0.4, scale$text_size / 12), "lines"),
      legend.key.width = ggplot2::unit(max(0.8, scale$text_size / 8), "lines")
    )
}

print_prob_heatmap_pdf <- function(model,
                                   data,
                                   plot.title,
                                   conformation,
                                   page_width = 11,
                                   page_height = 8.5) {
  if (!requireNamespace("grid", quietly = TRUE)) stop("`grid` package is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("`ggplot2` is required.")

  n_rows <- nrow(data)
  n_cols <- length(model$lev) + 1L
  scale <- prob_heatmap_scale_params(n_rows)
  p <- scale_prob_heatmap_plot(
    prob.heatmap(model, data, plot.title, conformation),
    n_rows = n_rows
  )
  grob <- ggplot2::ggplotGrob(p)

  cell_w <- 0.55
  cell_h <- cell_w * scale$plot_ratio
  plot_w_in <- min(page_width * 0.96, max(4, n_cols * cell_w + 1.0))
  plot_h_in <- min(page_height * 0.85, max(3, n_rows * cell_h + 2.0))

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(plot_w_in / page_width, "npc"),
    height = grid::unit(plot_h_in / page_height, "npc"),
    x = 0.5,
    y = 0.45,
    just = c("center", "center")
  ))
  grid::grid.draw(grob)
  grid::popViewport()
}

cv_classification <- function(data,
                              test.form,
                              ordinal = FALSE,
                              k = 5,
                              n.iter = 1,
                              seed = NULL,
                              verbose = FALSE,
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
      Test.data <- data[test_idx, , drop = FALSE]
      Train.data[[resp_var_clean]] <- factor(Train.data[[resp_var_clean]], levels = classes)

      if (length(unique(Train.data[[resp_var_clean]])) < 2) next

      if (!ordinal) {
        model <- nnet::multinom(test.form1, data = Train.data, maxit = 2000, trace = FALSE)
      } else {
        model <- fit_polr(formula = test.form, data = Train.data)
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
      names(row)[(ncol(row) - length(classes) + 1):ncol(row)] <- classes

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

# Deuteration (Kennard–Stone validation) --------------------------------------------

n_per_class <- NULL

data <- prepare_deuteration_training_data(training_csv)

train_idx <- build_kennard_stone_train_indices(
  data = data,
  class_col = "class",
  n_per_class = n_per_class
)
Train.data <- data[train_idx, , drop = FALSE]
Test.data <- data[-train_idx, , drop = FALSE]
External <- load_external_validation_aligned(external_csv, Train.data)

test.form <- "`class` ~ `X.2.3.` + `dip_y` + `diff.1.2` + `diff.2.3`"

log_msg("=== Setup (Kennard–Stone train pool for CV) ===")
cat(sprintf("Formula: %s\n", test.form))
cat(sprintf(
  "Samples: n_total=%d | KS train pool=%d | KS holdout=%d\n",
  nrow(data), nrow(Train.data), nrow(Test.data)
))
cat("Class counts (KS train pool):\n")
print(table(Train.data$class))
cat("Class counts (KS holdout):\n")
print(table(Test.data$class))
if (!is.null(External) && nrow(External) > 0) {
  cat(sprintf("External validation samples: n=%d\n", nrow(External)))
  cat("Class counts (external validation):\n")
  print(table(External$class))
} else {
  log_msg("External validation data not found or column layout mismatch; skipping external evaluation.")
}

CV7 <- cv_classification(
  Train.data,
  test.form,
  ordinal = TRUE,
  k = 7,
  n.iter = 200,
  seed = 10
)

LOO <- cv_classification(
  Train.data,
  test.form,
  ordinal = TRUE,
  k = "loo",
  n.iter = 200,
  seed = 10
)

log_msg("=== Cross-validation performance (on Kennard–Stone training pool) ===")
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

report_dir <- file.path(dirname(dirname(script_dir)), "Model Reports")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)
report_file <- file.path(report_dir, "CrossValidation.ModelReport.KennardStone.pdf")

log_msg(sprintf("Writing cross-validation model report PDF: %s", report_file))

pdf(report_file, width = 11, height = 8.5)

acc_cv_pct <- CV7$results_table$accuracy * 100
acc_loo_pct <- LOO$results_table$accuracy * 100

par(mar = c(5, 4, 4, 2) + 0.1)
boxplot(
  list("CV (k=7)" = acc_cv_pct, "LOO" = acc_loo_pct),
  main = "Cross-validation accuracy (Kennard–Stone train pool)",
  ylab = "Accuracy (%)",
  col = c("#4C78A8", "#F58518")
)

hist(acc_cv_pct, breaks = 20, col = "#4C78A8AA", border = "#4C78A8",
     main = "CV (k=7) accuracy histogram", xlab = "Accuracy (%)")
hist(acc_loo_pct, breaks = 20, col = "#F58518AA", border = "#F58518",
     main = "LOO accuracy histogram", xlab = "Accuracy (%)")

final_model <- fit_polr(formula = test.form, data = Train.data)

log_msg("=== Final fitted model (Kennard–Stone train pool) ===")
print(summary(final_model))
cat(sprintf("Accuracy on KS train pool: %.2f%%\n", pred_accuracy(final_model, Train.data) * 100))
if (nrow(Test.data) > 0) {
  cat(sprintf("Accuracy on KS holdout: %.2f%%\n", pred_accuracy(final_model, Test.data) * 100))
}
mi_ex <- if (!is.null(External) && nrow(External) > 0) {
  cat(sprintf("Accuracy on external validation: %.2f%%\n", pred_accuracy(final_model, External) * 100))
  mod.info(final_model, External, FALSE, FALSE)
} else {
  NULL
}

mi_tr <- mod.info(final_model, Train.data, FALSE, FALSE)
mi_te <- if (nrow(Test.data) > 0) mod.info(final_model, Test.data, FALSE, FALSE) else NULL

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

grid_text_table_page <- function(title, table_lines) {
  if (!requireNamespace("grid", quietly = TRUE)) stop("`grid` package is required.")
  grid::grid.newpage()
  grid::grid.text(title, x = 0.5, y = 0.98, just = "top",
                  gp = grid::gpar(fontface = "bold", cex = 1.0))
  grid::grid.text(table_lines, x = 0.02, y = 0.94, just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "mono", cex = 0.55),
                  default.units = "npc")
}

cm_train <- ct_plot(
  mi_tr$class.table,
  plot.title = "Train Set (Kennard–Stone pool; model report)",
  conformation = "3. 3rd Place"
)
grid_side_by_side_confusion(
  confusion_plot = cm_train$plot,
  confusion_table = mi_tr$class.table,
  left_title = "Train confusion matrix",
  right_title = "Train confusion table"
)

print_prob_heatmap_pdf(
  final_model,
  Train.data,
  plot.title = "Train Set (Kennard–Stone pool; model report)",
  conformation = "3. 3rd Place"
)

if (!is.null(mi_te) && nrow(Test.data) > 0) {
  cm_test <- ct_plot(
    mi_te$class.table,
    plot.title = "Holdout (Kennard–Stone; model report)",
    conformation = "3. 3rd Place"
  )
  grid_side_by_side_confusion(
    confusion_plot = cm_test$plot,
    confusion_table = mi_te$class.table,
    left_title = "Holdout confusion matrix",
    right_title = "Holdout confusion table"
  )
  print_prob_heatmap_pdf(
    final_model,
    Test.data,
    plot.title = "Holdout (Kennard–Stone; model report)",
    conformation = "3. 3rd Place"
  )
}

if (!is.null(mi_ex) && !is.null(External) && nrow(External) > 0) {
  ext_pred_tbl <- build_prediction_table(final_model, External)
  ext_pred_display <- data.frame(
    substrate = rownames(External),
    ext_pred_tbl,
    check.names = FALSE
  )
  rownames(ext_pred_display) <- NULL
  num_cols <- vapply(ext_pred_display, is.numeric, logical(1))
  ext_pred_display[num_cols] <- lapply(ext_pred_display[num_cols], round, 2)
  ext_pred_lines <- capture.output(print(ext_pred_display, row.names = FALSE))
  grid_text_table_page(
    title = sprintf(
      "External validation — class probabilities and predictions (accuracy = %.2f%%)",
      pred_accuracy(final_model, External) * 100
    ),
    table_lines = paste(ext_pred_lines, collapse = "\n")
  )

  cm_ext <- ct_plot(
    mi_ex$class.table,
    plot.title = "External Validation (model report)",
    conformation = "1. 1st Place"
  )
  grid_side_by_side_confusion(
    confusion_plot = cm_ext$plot,
    confusion_table = mi_ex$class.table,
    left_title = "External confusion matrix",
    right_title = "External confusion table"
  )
  print_prob_heatmap_pdf(
    final_model,
    External,
    plot.title = "External Validation (model report)",
    conformation = "1. 1st Place"
  )
}

dev.off()
log_msg("Cross-validation model report PDF written successfully.")
log_msg("Finished Cross validation_Kennard_Stone.R.")
