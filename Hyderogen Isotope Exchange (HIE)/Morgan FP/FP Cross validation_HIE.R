library("rxn.cond.class")

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
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

scaled_prob_heatmap <- function(model, data, plot.title, conformation) {
  scale_prob_heatmap_plot(
    prob.heatmap(model, data, plot.title, conformation),
    n_rows = nrow(data)
  )
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

rank_conformation <- function(rank) {
  rank <- as.integer(rank)
  suffix <- if (rank %% 100L %in% 11:13) {
    "th"
  } else {
    c("th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th")[rank %% 10L + 1L]
  }
  sprintf("%d. %d%s Place", rank, rank, suffix)
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

# Cross validation code ---------------------------------------------------

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

# Hydrogen Isotope Exchange — Morgan FP (redundant bits removed)

.morgan_dir <- .morgan_features_dir()
source(file.path(.morgan_dir, "morgan_feature_config.R"), local = TRUE)
.morgan_load <- load_morgan_feature_table(.morgan_dir)
.morgan_spec <- .morgan_load$spec
training_csv <- .morgan_load$training_csv
feat_cols <- .morgan_load$feature_cols
data50 <- .morgan_load$data

## Sampler data (feature space for similarity; same row order as data50)
sampl.dat <- data50[, c(feat_cols, "class", "flag"), drop = FALSE]
sampl.dat$class <- factor(sampl.dat$class)

one <- simi.sampler(sampl.dat, 1)
two <- simi.sampler(sampl.dat, 2)
three <- simi.sampler(sampl.dat, 3)
four <- simi.sampler(sampl.dat, 4)

two_three <- simi.sampler(sampl.dat, 2, 3)
one_three <- simi.sampler(sampl.dat, 1, 3)
four_three <- simi.sampler(sampl.dat, 4, 3)

similarties <- c(
  three,
  union(two, two_three),
  union(one, one_three),
  union(four, four_three)
)

Train.data <- data50[similarties, , drop = FALSE]
Test.data <- data50[-similarties, , drop = FALSE]

# -------------------------------------------------------------------------
# User-defined model (edit after FP Model Construction_HIE.R)
# -------------------------------------------------------------------------

chosen_rank <- 10L
test.form <- "`class` ~ `MFP_448` + `MFP_694` + `MFP_807` + `MFP_1160`"
model_conformation <- rank_conformation(chosen_rank)

log_msg(sprintf("=== User-defined model (%s) ===", .morgan_spec$label))
cat(sprintf("Chosen rank: %d (%s)\n", chosen_rank, model_conformation))
cat(sprintf("Chosen formula: %s\n", test.form))

CV7 <- cv_classification(
  Train.data,
  test.form,
  ordinal = FALSE,
  k = 7,
  n.iter = 200,
  seed = 100
)
message("7-fold CV mean accuracy (%): ", round(CV7$overall_mean_accuracy * 100, 2))

LOO <- cv_classification(
  Train.data,
  test.form,
  ordinal = FALSE,
  k = "loo",
  n.iter = 1,
  seed = 100
)
message("LOO mean accuracy (%): ", round(LOO$overall_mean_accuracy * 100, 2))

# Fit once on training set and report performance on train and held-out test
final_model <- nnet::multinom(test.form, data = Train.data, maxit = 2000, trace = FALSE)

train_pred <- predict(final_model, newdata = Train.data, type = "class")
test_pred <- predict(final_model, newdata = Test.data, type = "class")

train_acc <- mean(as.character(train_pred) == as.character(Train.data$class))
test_acc <- mean(as.character(test_pred) == as.character(Test.data$class))

train_mcfadden_r2 <- calc_mcfadden_r2(final_model, Train.data, "class")
test_mcfadden_r2 <- calc_mcfadden_r2(final_model, Test.data, "class")

message("Fixed formula: ", deparse(test.form))
message("Train accuracy (%): ", round(train_acc * 100, 2))
message("Test accuracy (%): ", round(test_acc * 100, 2))
message("Train McFadden R2: ", round(train_mcfadden_r2, 3))
message("Test McFadden R2: ", round(test_mcfadden_r2, 3))
message("Train confusion matrix:")
print(table(actual = Train.data$class, predicted = train_pred))
message("Test confusion matrix:")
print(table(actual = Test.data$class, predicted = test_pred))

# -------------------------------
# Model report (PDF) generation
# -------------------------------
report_dir <- .morgan_dir
report_file <- file.path(report_dir, "CrossValidation.ModelReport.HIE.MorganFP.pdf")

log_msg(sprintf("Writing cross-validation model report PDF: %s", report_file))

pdf(report_file, width = 11, height = 8.5)

# Accuracy distributions from CV / LOO
acc_cv_pct <- CV7$results_table$accuracy * 100
acc_loo_pct <- LOO$results_table$accuracy * 100

par(mar = c(5, 4, 4, 2) + 0.1)
boxplot(
  list("CV (k=7)" = acc_cv_pct, "LOO" = acc_loo_pct),
  main = paste0("Cross-validation accuracy distributions (HIE ", .morgan_spec$label, ")"),
  ylab = "Accuracy (%)",
  col = c("#4C78A8", "#F58518")
)

hist(acc_cv_pct, breaks = 20, col = "#4C78A8AA", border = "#4C78A8",
     main = "CV (k=7) accuracy histogram", xlab = "Accuracy (%)")
hist(acc_loo_pct, breaks = 20, col = "#F58518AA", border = "#F58518",
     main = "LOO accuracy histogram", xlab = "Accuracy (%)")

# Text summary page (formula + key metrics)
plot.new()
summary_lines <- c(
  paste0("Hydrogen Isotope Exchange - ", .morgan_spec$label),
  paste0("Chosen rank: ", chosen_rank, " (", model_conformation, ")"),
  paste0("Chosen formula: ", test.form),
  paste0("Train accuracy (%): ", round(train_acc * 100, 2)),
  paste0("Test accuracy (%): ", round(test_acc * 100, 2)),
  paste0("Train McFadden R2: ", round(train_mcfadden_r2, 3)),
  paste0("Q^2 (McFadden on test): ", round(test_mcfadden_r2, 3)),
  paste0("CV (k=7) mean accuracy (%): ", round(CV7$overall_mean_accuracy * 100, 2)),
  paste0("LOO mean accuracy (%): ", round(LOO$overall_mean_accuracy * 100, 2))
)
text(0.02, 0.98, paste(summary_lines, collapse = "\n"), adj = c(0, 1), cex = 1.0)

grid_side_by_side_plots <- function(left_plot,
                                    right_plot,
                                    left_title,
                                    right_title,
                                    right_n_rows = NULL,
                                    right_n_cols = NULL,
                                    page_width = 11,
                                    page_height = 8.5) {
  if (!requireNamespace("grid", quietly = TRUE)) stop("`grid` package is required.")

  as_grob <- function(p) {
    if (inherits(p, "ggplot")) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("`ggplot2` is required to render ggplot objects.")
      }
      return(ggplot2::ggplotGrob(p))
    }
    p
  }

  left_grob <- as_grob(left_plot)
  right_grob <- as_grob(right_plot)

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))

  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(left_grob)
  grid::grid.text(left_title, x = 0.5, y = 0.98, just = "top",
                  gp = grid::gpar(fontface = "bold", cex = 1.0))
  grid::popViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.text(right_title, x = 0.5, y = 0.98, just = "top",
                  gp = grid::gpar(fontface = "bold", cex = 1.0))
  if (!is.null(right_n_rows)) {
    scale <- prob_heatmap_scale_params(right_n_rows)
    n_cols <- if (is.null(right_n_cols)) 4L else right_n_cols
    cell_w <- 0.55
    cell_h <- cell_w * scale$plot_ratio
    panel_w_in <- page_width / 2
    plot_w_in <- min(panel_w_in * 0.92, max(2.5, n_cols * cell_w + 0.8))
    plot_h_in <- min(page_height * 0.85, max(2.5, right_n_rows * cell_h + 1.5))
    grid::pushViewport(grid::viewport(
      width = grid::unit(plot_w_in / panel_w_in, "npc"),
      height = grid::unit(plot_h_in / page_height, "npc"),
      x = 0.5,
      y = 0.45,
      just = c("center", "center")
    ))
    grid::grid.draw(right_grob)
    grid::popViewport()
  } else {
    grid::grid.draw(right_grob)
  }
  grid::popViewport()

  grid::popViewport()
}

model.info.train <- mod.info(final_model, Train.data, FALSE, FALSE)
cm_train <- ct_plot(
  model.info.train$class.table,
  plot.title = "Train Set (for model report)",
  conformation = model_conformation
)
p_train_hm <- scaled_prob_heatmap(
  final_model,
  Train.data,
  plot.title = "Train Set (for model report)",
  conformation = model_conformation
)
grid_side_by_side_plots(
  left_plot = cm_train$plot,
  right_plot = p_train_hm,
  left_title = "Train confusion matrix",
  right_title = "Train probability heatmap",
  right_n_rows = nrow(Train.data),
  right_n_cols = length(final_model$lev) + 1L
)

if (nrow(Test.data) > 0) {
  model.info.test <- mod.info(final_model, Test.data, FALSE, FALSE)
  cm_test <- ct_plot(
    model.info.test$class.table,
    plot.title = "Test Set (for model report)",
    conformation = model_conformation
  )
  p_test_hm <- scaled_prob_heatmap(
    final_model,
    Test.data,
    plot.title = "Test Set (for model report)",
    conformation = model_conformation
  )
  grid_side_by_side_plots(
    left_plot = cm_test$plot,
    right_plot = p_test_hm,
    left_title = "Test confusion matrix",
    right_title = "Test probability heatmap",
    right_n_rows = nrow(Test.data),
    right_n_cols = length(final_model$lev) + 1L
  )
}

dev.off()
log_msg("Cross-validation model report PDF written successfully.")
