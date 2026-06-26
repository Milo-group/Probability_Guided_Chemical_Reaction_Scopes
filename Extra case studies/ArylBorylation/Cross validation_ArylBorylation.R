library('rxn.cond.class')

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) > 0 && nzchar(file_path[1])) {
    return(dirname(normalizePath(file_path[1], winslash = "/", mustWork = FALSE)))
  }

  src_path <- NULL
  src_path <- tryCatch(
    normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE),
    error = function(e) NULL
  )
  if (!is.null(src_path) && nzchar(src_path)) {
    return(dirname(src_path))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

make_sample_ids <- function(n) {
  # For this case (22 samples), IDs will be A..V.
  # If n > 26, continue with S27, S28, ... to keep uniqueness.
  if (n <= 26) {
    return(LETTERS[seq_len(n)])
  }
  c(LETTERS, paste0("S", seq(27, n)))
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
  
  # Accept either a formula object or a formula string.
  # If the input string incorrectly uses multiple "~", keep the first one
  # and convert the rest to "+" so R can parse it as a valid model formula.
  if (is.character(test.form)) {
    form_txt <- trimws(test.form)
    if (length(gregexpr("~", form_txt, fixed = TRUE)[[1]]) > 1) {
      parts <- strsplit(form_txt, "~", fixed = TRUE)[[1]]
      lhs <- trimws(parts[1])
      rhs_terms <- trimws(parts[-1])
      form_txt <- paste(lhs, "~", paste(rhs_terms, collapse = " + "))
      message("Normalized formula: ", form_txt)
    }
    test.form1 <- as.formula(form_txt)
  } else {
    test.form1 <- test.form
  }

  # Extract response variable robustly from formula
  resp_var_clean <- all.vars(test.form1)[1]
  
  # Check response column exists and factorize
  if (!resp_var_clean %in% names(data)) stop("response variable not found in data")
  if (!is.factor(data[[resp_var_clean]])) {
    warning("Response variable coerced to factor for stratified sampling.")
    data[[resp_var_clean]] <- factor(data[[resp_var_clean]])
  }
  
  n <- nrow(data)
  classes <- levels(data[[resp_var_clean]])
  sample_ids <- rownames(data)
  if (is.null(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    sample_ids <- make_sample_ids(n)
  }
  
  # Handle "LOO"
  loo <- FALSE
  if (is.character(k) && tolower(k) == "loo") {
    k <- n
    loo <- TRUE
  }
  
  # Prevent repeats in LOO
  if (k == n && n.iter > 1) {
    if (verbose) message("LOO detected. Repeats ignored; performing LOO once.")
    n.iter <- 1
  }
  
  results_table <- data.frame()
  
  # Progress bar
  total_folds <- n.iter * k
  pb <- utils::txtProgressBar(min = 0, max = total_folds, style = 3)
  progress_counter <- 0
  
  for (iter in seq_len(n.iter)) {
    if (!is.null(seed)) set.seed(seed + iter)
    
    if (k == n) {
      folds <- as.list(seq_len(n))  # LOO
    } else {
      folds <- caret::createFolds(data[[resp_var_clean]], k = k, list = TRUE, returnTrain = FALSE)
    }
    
    for (f in seq_along(folds)) {
      test_idx <- folds[[f]]
      Train.data <- data[-test_idx, , drop = FALSE]
      Test.data  <- data[test_idx, , drop = FALSE]
      
      # Ensure all factor levels are retained in training fold
      Train.data[[resp_var_clean]] <- factor(Train.data[[resp_var_clean]], levels = classes)
      
      # Skip if <2 classes in training
      if (length(unique(Train.data[[resp_var_clean]])) < 2) next
      
      # Fit model (nnet::multinom uses a random initialization of weights, so LOO my alter)
      if (!ordinal) {
        model <- nnet::multinom(test.form1, data = Train.data, maxit = 2000, trace = FALSE)
      } else {
        model <- fit_polr(formula = test.form1, data = Train.data)
      }
      
      # Predict class
      pred_class <- predict(model, newdata = Test.data, type = "class")
      acc <- mean(as.character(pred_class) == as.character(Test.data[[resp_var_clean]]))
      
      # Class counts
      class_counts <- table(factor(Test.data[[resp_var_clean]], levels = classes))
      
      # Build row
      row <- data.frame(
        iteration = iter,
        fold = f,
        left_out_samples = paste(sample_ids[test_idx], collapse = ","),
        accuracy = acc
      )
      
      # Add predicted class only for LOO
      if (loo) row$predicted_class <- as.character(pred_class)
      
      # Add class counts
      row <- cbind(row, as.list(as.numeric(class_counts)))
      names(row)[(ncol(row)-length(classes)+1):ncol(row)] <- classes
      
      results_table <- rbind(results_table, row)
      
      # Update progress bar (always shown)
      progress_counter <- progress_counter + 1
      utils::setTxtProgressBar(pb, progress_counter)
      
      # Verbose message (optional)
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

# Aryl Borylation ---------------------------------------------------------
# Load data
.aryl_dir <- get_script_dir()
ArylBorylation <- read.csv(
  file.path(.aryl_dir, 'Data', 'Stevens_data_organized_for_classification_fixed.csv'),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(ArylBorylation) <- ArylBorylation[['Product Name']]

# Remove identifier / metadata columns by name
drop_cols <- c(
  'Electrophile', 'Electrophile_inchi',
  'Ligand', 'Ligand_inchi',
  'Product_inchi', 'Product Name old', 'Product Name',
  'Solvent_inchi', 'Yield'
)
ArylBorylation <- ArylBorylation[, !(names(ArylBorylation) %in% drop_cols)]

# Remove zero-variance numeric columns
is_zero_var <- vapply(
  ArylBorylation,
  function(x) is.numeric(x) && (length(unique(x)) == 1),
  logical(1)
)
ArylBorylation <- ArylBorylation[, !is_zero_var]

# Ensure numeric predictors and factor response
predictor_cols <- setdiff(names(ArylBorylation), 'class')
numeric_predictors <- lapply(
  ArylBorylation[predictor_cols],
  function(x) suppressWarnings(as.numeric(as.character(x)))
)
valid_predictors <- vapply(
  predictor_cols,
  function(col_name) {
    original <- ArylBorylation[[col_name]]
    converted <- numeric_predictors[[col_name]]
    all(is.na(converted) == is.na(original))
  },
  logical(1)
)

if (!all(valid_predictors)) {
  message('Dropping non-numeric predictors: ',
          paste(predictor_cols[!valid_predictors], collapse = ', '))
}

predictor_cols <- predictor_cols[valid_predictors]
ArylBorylation[predictor_cols] <- numeric_predictors[predictor_cols]
ArylBorylation$class <- suppressWarnings(as.numeric(as.character(ArylBorylation$class)))

# Keep only complete finite rows before similarity sampling
keep_rows <- complete.cases(ArylBorylation[, c(predictor_cols, 'class'), drop = FALSE]) &
  apply(ArylBorylation[, predictor_cols, drop = FALSE], 1, function(x) all(is.finite(x)))
ArylBorylation <- ArylBorylation[keep_rows, , drop = FALSE]
row.names(ArylBorylation) <- make_sample_ids(nrow(ArylBorylation))

ArylBorylation$class <- as.factor(ArylBorylation$class)

# Add a 'flag' column to sequentially number the rows
ArylBorylation <- plyr::mutate(ArylBorylation, flag = seq(1,nrow(ArylBorylation)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(ArylBorylation, 1, sample.size = 5)  # Sample from class 1
two <- simi.sampler(ArylBorylation, 2, sample.size = 5)  # Sample from class 2
three <- simi.sampler(ArylBorylation, 3, sample.size = 5)  # Sample from class 3

# # Sample class 1 molecules similar to class 3
# one_three <- simi.sampler(ArylBorylation, 1, 3)  

# # Sample class 2 molecules similar to class 3
# two_three <- simi.sampler(ArylBorylation, 2, 3)  

# Combine the similarities from the various classes into one vector
# similarties <- c(union(one, one_three), 
#                  union(two, two_three), 
#                  three)
                 
similarties <- c(one,two,three)
# Define train and test data from the samples taken from groups 1 and 2
Train.data <- ArylBorylation[similarties,]
Test.data <- ArylBorylation[-similarties,]

test.form <- "`class` ~ `ArylHalide_XY_Length_Min` + `ArylHalide_XY_Lowdin_Max` + `ArylHalide_XY_Mulliken_Bond_Max`"
test.formula <- as.formula(test.form)
# `ArylHalide_C_Chem Shift`| ArylHalide_C_Mulliken| `ArylHalide_XY_Exposed Area_Avg`| ArylHalide_XY_Mulliken_Max
# Train/test performance after similarity sampling
model_after_simi <- nnet::multinom(test.formula,
                                   data = Train.data,
                                   maxit = 2000,
                                   trace = FALSE)
null_model_after_simi <- nnet::multinom(class ~ 1,
                                        data = Train.data,
                                        maxit = 2000,
                                        trace = FALSE)

class_levels <- levels(Train.data$class)
train_pred <- predict(model_after_simi, newdata = Train.data, type = "class")
test_pred <- predict(model_after_simi, newdata = Test.data, type = "class")

train_acc <- mean(as.character(train_pred) == as.character(Train.data$class))
test_acc <- mean(as.character(test_pred) == as.character(Test.data$class))
loglik_model <- as.numeric(stats::logLik(model_after_simi))
loglik_null <- as.numeric(stats::logLik(null_model_after_simi))
mcfadden_r2 <- 1 - (loglik_model / loglik_null)

train_cm <- table(
  Observed = factor(Train.data$class, levels = class_levels),
  Predicted = factor(train_pred, levels = class_levels)
)
test_cm <- table(
  Observed = factor(Test.data$class, levels = class_levels),
  Predicted = factor(test_pred, levels = class_levels)
)

cat("\nAfter similarity sampling - overall accuracy (%):\n")
cat("Train:", round(train_acc * 100, 2), "\n")
cat("Test :", round(test_acc * 100, 2), "\n")
cat("McFadden R2:", round(mcfadden_r2, 4), "\n")

cat("\nAfter similarity sampling - confusion matrix (Train):\n")
print(train_cm)

cat("\nAfter similarity sampling - confusion matrix (Test):\n")
print(test_cm)

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

cat("\n7-fold CV mean accuracy (%):\n")
print(round(CV7$overall_mean_accuracy * 100, 2))
cat("\n7-fold CV results (first 10 rows):\n")
print(utils::head(CV7$results_table, 10))

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)

cat("\nLOO CV mean accuracy (%):\n")
print(round(LOO$overall_mean_accuracy * 100, 2))
cat("\nLOO CV results (first 10 rows):\n")
print(utils::head(LOO$results_table, 10))

# ---------------------------------
# Model report (PDF) generation
# ---------------------------------
report_dir <- file.path(dirname(.aryl_dir), "Model Reports")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)
report_file <- file.path(report_dir, "CrossValidation.ModelReport.ArylBorylation.pdf")

log_msg(sprintf("Writing cross-validation model report PDF: %s", report_file))

pdf(report_file, width = 11, height = 8.5)

# Accuracy distributions from CV / LOO
acc_cv_pct <- CV7$results_table$accuracy * 100
acc_loo_pct <- LOO$results_table$accuracy * 100

par(mar = c(5, 4, 4, 2) + 0.1)
boxplot(
  list("CV (k=7)" = acc_cv_pct, "LOO" = acc_loo_pct),
  main = "Cross-validation accuracy distributions (Aryl Borylation)",
  ylab = "Accuracy (%)",
  col = c("#4C78A8", "#F58518")
)

hist(
  acc_cv_pct, breaks = 20, col = "#4C78A8AA", border = "#4C78A8",
  main = "CV (k=7) accuracy histogram", xlab = "Accuracy (%)"
)
hist(
  acc_loo_pct, breaks = 20, col = "#F58518AA", border = "#F58518",
  main = "LOO accuracy histogram", xlab = "Accuracy (%)"
)

# Text summary page (formula + key metrics)
plot.new()
summary_lines <- c(
  "Aryl Borylation",
  paste0("Formula: ", test.form),
  paste0("Train accuracy (%): ", round(train_acc * 100, 2)),
  paste0("Test accuracy (%): ", round(test_acc * 100, 2)),
  paste0("McFadden R2: ", round(mcfadden_r2, 3)),
  paste0("CV (k=7) mean accuracy (%): ", round(CV7$overall_mean_accuracy * 100, 2)),
  paste0("LOO mean accuracy (%): ", round(LOO$overall_mean_accuracy * 100, 2))
)
text(0.02, 0.98, paste(summary_lines, collapse = "\n"), adj = c(0, 1), cex = 1.0)

grid_side_by_side_plots <- function(left_plot, right_plot, left_title, right_title) {
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
  grid::grid.draw(right_grob)
  grid::grid.text(right_title, x = 0.5, y = 0.98, just = "top",
                  gp = grid::gpar(fontface = "bold", cex = 1.0))
  grid::popViewport()

  grid::popViewport()
}

model.info.train <- mod.info(model_after_simi, Train.data, FALSE, FALSE)
cm_train <- ct_plot(
  model.info.train$class.table,
  plot.title = "Train Set (for model report)",
  conformation = "1. 1st Place"
)
p_train_hm <- prob.heatmap(
  model_after_simi,
  Train.data,
  plot.title = "Train Set (for model report)",
  conformation = "1. 1st Place"
)
grid_side_by_side_plots(
  left_plot = cm_train$plot,
  right_plot = p_train_hm,
  left_title = "Train confusion matrix",
  right_title = "Train probability heatmap"
)

if (nrow(Test.data) > 0) {
  model.info.test <- mod.info(model_after_simi, Test.data, FALSE, FALSE)
  cm_test <- ct_plot(
    model.info.test$class.table,
    plot.title = "Test Set (for model report)",
    conformation = "1. 1st Place"
  )
  p_test_hm <- prob.heatmap(
    model_after_simi,
    Test.data,
    plot.title = "Test Set (for model report)",
    conformation = "1. 1st Place"
  )
  grid_side_by_side_plots(
    left_plot = cm_test$plot,
    right_plot = p_test_hm,
    left_title = "Test confusion matrix",
    right_title = "Test probability heatmap"
  )
}

dev.off()
log_msg("Cross-validation model report PDF written successfully.")
