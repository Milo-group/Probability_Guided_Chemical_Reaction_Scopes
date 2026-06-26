# Getting Started Example ג€” cross-validation
# Mirrors Aldehyde-Deuteration/Data/Cross validation_Deuteration.R (lighter settings for speed).

library("rxn.cond.class")

get_script_dir <- function() {
  path <- tryCatch(
    normalizePath(sys.frame(1)$ofile, winslash = "/"),
    error = function(e) NA_character_
  )
  if (!is.na(path)) return(dirname(path))
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  if (nzchar(file_arg)) return(dirname(normalizePath(file_arg, winslash = "/")))
  "."
}

script_dir <- get_script_dir()
training_csv <- file.path(script_dir, "Training_Data.csv")

cat("=== Step 1: Load data and create train / test split ===\n")
data <- data.frame(data.table::fread(training_csv), check.names = FALSE)
row.names(data) <- data[, 2]
data$class <- as.factor(data$class)
data <- data[, -c(1:2)]

one <- simi.sampler(data, 1, sample.size = 6)
two <- simi.sampler(data, 2, sample.size = 6)
three <- simi.sampler(data, 3, sample.size = 6)
one_three <- simi.sampler(data, 1, 3, sample.size = 3)
two_three <- simi.sampler(data, 2, 3, sample.size = 3)
similarties <- unique(c(union(one, one_three), union(two, two_three), three))

Train.data <- data[similarties, ]
Test.data <- data[-similarties, ]

cat("\n=== Step 2: Quick model search ===\n")
models.ordinal <- sub_model_log(data = Train.data, min = 2, max = 2, ordinal = TRUE)
test.form <- models.ordinal[1, 1][[1]]
cat("Using formula:", gsub("`", "", test.form), "\n")

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
    data[[resp_var_clean]] <- factor(data[[resp_var_clean]])
  }

  n <- nrow(data)
  classes <- levels(data[[resp_var_clean]])

  loo <- FALSE
  if (is.character(k) && tolower(k) == "loo") {
    k <- n
    loo <- TRUE
  }
  if (k == n && n.iter > 1) n.iter <- 1

  results_table <- data.frame()
  total_folds <- n.iter * if (loo) n else k
  pb <- utils::txtProgressBar(min = 0, max = total_folds, style = 3)
  progress_counter <- 0

  for (iter in seq_len(n.iter)) {
    if (!is.null(seed)) set.seed(seed + iter)
    folds <- if (k == n) as.list(seq_len(n)) else caret::createFolds(data[[resp_var_clean]], k = k, list = TRUE, returnTrain = FALSE)

    for (f in seq_along(folds)) {
      test_idx <- folds[[f]]
      Train.fold <- data[-test_idx, , drop = FALSE]
      Test.fold <- data[test_idx, , drop = FALSE]
      Train.fold[[resp_var_clean]] <- factor(Train.fold[[resp_var_clean]], levels = classes)
      if (length(unique(Train.fold[[resp_var_clean]])) < 2) next

      model <- if (!ordinal) {
        nnet::multinom(test.form1, data = Train.fold, maxit = 2000, trace = FALSE)
      } else {
        fit_polr(formula = test.form, data = Train.fold)
      }

      pred_class <- predict(model, newdata = Test.fold, type = "class")
      acc <- mean(as.character(pred_class) == as.character(Test.fold[[resp_var_clean]]))
      class_counts <- table(factor(Test.fold[[resp_var_clean]], levels = classes))
      row <- data.frame(iteration = iter, fold = f, accuracy = acc)
      row <- cbind(row, as.list(as.numeric(class_counts)))
      names(row)[(ncol(row) - length(classes) + 1):ncol(row)] <- classes
      results_table <- rbind(results_table, row)

      progress_counter <- progress_counter + 1
      utils::setTxtProgressBar(pb, progress_counter)
    }
  }
  close(pb)

  list(
    results_table = results_table,
    overall_mean_accuracy = mean(results_table$accuracy, na.rm = TRUE),
    k = if (k == n) "LOO" else k,
    n.iter = n.iter
  )
}

cat("\n=== Step 3: Cross-validation (5-fold, 10 repeats) ===\n")
CV5 <- cv_classification(Train.data, test.form, ordinal = TRUE, k = 5, n.iter = 10, seed = 10)
cat(sprintf("\nMean CV accuracy: %.1f%%\n", CV5$overall_mean_accuracy * 100))

cat("\n=== Step 4: Leave-one-out cross-validation ===\n")
LOO <- cv_classification(Train.data, test.form, ordinal = TRUE, k = "loo", n.iter = 1, seed = 10)
cat(sprintf("LOO accuracy: %.1f%%\n", LOO$overall_mean_accuracy * 100))

cat("\n=== Step 5: Final model on full training subset ===\n")
final_model <- fit_polr(formula = test.form, data = Train.data)
mi_train <- mod.info(final_model, Train.data, TRUE, FALSE)
mi_test <- mod.info(final_model, Test.data, TRUE, FALSE)
cat(sprintf("Train accuracy: %.1f%% | Test accuracy: %.1f%%\n", mi_train$accuracy, mi_test$accuracy))

cat("\nDone. For a full PDF report workflow, see Aldehyde-Deuteration/Data/Cross validation_Deuteration.R\n")
