library('rxn.cond.class')


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
  
  # Keep original string for fit_polr
  # Convert to formula only for extracting response
  test.form1 <- as.formula(test.form)
  
  # Extract response variable safely from formula object
  resp_var <- deparse(test.form1[[2]])          # turns objecct to string          
  resp_var_clean <- gsub("`", "", resp_var)     # clean name for data indexing
  
  # Check response column exists and factorize
  if (!resp_var_clean %in% names(data)) stop("response variable not found in data")
  if (!is.factor(data[[resp_var_clean]])) {
    warning("Response variable coerced to factor for stratified sampling.")
    data[[resp_var_clean]] <- factor(data[[resp_var_clean]])
  }
  
  n <- nrow(data)
  classes <- levels(data[[resp_var_clean]])
  
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
        # Ordinal: use the original string formula with backticks
        model <- fit_polr(formula = test.form, data = Train.data)
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
        left_out_samples = paste(test_idx, collapse = ","),
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

# Deuteration -------------------------------------------------------------

# Load data from a CSV file into a data frame
data <- data.frame(data.table::fread(training_csv), check.names = F)

# Clean and organize data
row.names(data) <- data[,2]  # Set the second column as row names
data$class <- as.factor(data$class)  # Convert the 'class' column to a factor
data <- data[,-c(1:2)]  # Remove the first and second columns (name and tag)

# Similarity-based sampling for each class
one <- simi.sampler(data, 1)  # Sample from class 1
two <- simi.sampler(data, 2)  # Sample from class 2
three <- simi.sampler(data, 3)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(data, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(data, 2, 3)  

# Combine the similarities from the various classes into one vector
similarties <- c(union(one, one_three), 
                 union(two, two_three), 
                 three)

# Define the training set from the subset of similarities
Train.data <- data[similarties, ]

# Define the test set from the remaining rows not in the similarities set
Test.data <- data[-similarties, ]

test.form <- "`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`"

# 7-fold cross-validation
CV7 <- cv_classification(Train.data, test.form, ordinal = TRUE,k = 7, 
                         n.iter = 200, seed = 10)
round(CV7$overall_mean_accuracy*100,2)
"[1] 79.22"
# Leave-one-out cross-validation
LOO <- cv_classification(Train.data, test.form, ordinal = TRUE,k = "loo", 
                         n.iter = 200, seed = 10)
round(LOO$overall_mean_accuracy*100,2)
"[1] 76.19"


# -------------------------------
# Model report (PDF) generation
# -------------------------------
report_dir <- file.path(dirname(script_dir), "Model Reports")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)
report_file <- file.path(report_dir, "CrossValidation.ModelReport.Deuteration.trial.pdf")

log_msg(sprintf("Writing cross-validation model report PDF: %s", report_file))

pdf(report_file, width = 11, height = 8.5)

# Accuracy distributions from CV / LOO
acc_cv_pct <- CV7$results_table$accuracy * 100
acc_loo_pct <- LOO$results_table$accuracy * 100

par(mar = c(5, 4, 4, 2) + 0.1)
boxplot(
  list("CV (k=7)" = acc_cv_pct, "LOO" = acc_loo_pct),
  main = "Cross-validation accuracy distributions",
  ylab = "Accuracy (%)",
  col = c("#4C78A8", "#F58518")
)

hist(acc_cv_pct, breaks = 20, col = "#4C78A8AA", border = "#4C78A8",
     main = "CV (k=7) accuracy histogram", xlab = "Accuracy (%)")
hist(acc_loo_pct, breaks = 20, col = "#F58518AA", border = "#F58518",
     main = "LOO accuracy histogram", xlab = "Accuracy (%)")


# Fit one final model on the full training subset, then create the usual reports.
final_model <- fit_polr(formula = test.form, data = Train.data)

grid_side_by_side_confusion <- function(confusion_plot, confusion_table, left_title, right_title) {
  # Draw the confusion plot (left) next to a text table (right) using base `grid`.
  # This avoids relying on external packages like `gridExtra`.
  if (!requireNamespace("grid", quietly = TRUE)) stop("`grid` package is required.")
  
  # Convert ggplot to grob when needed; otherwise assume it's already a grob.
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
  
  # Left: confusion matrix graphic
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(confusion_grob)
  grid::grid.text(left_title, x = 0.5, y = 0.98, just = "top",
                  gp = grid::gpar(fontface = "bold", cex = 1.0))
  grid::popViewport()
  
  # Right: numeric confusion matrix table
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.text(right_title, x = 0.5, y = 0.98, just = "top",
                  gp = grid::gpar(fontface = "bold", cex = 1.0))
  grid::grid.text(tab_label, x = 0.02, y = 0.96, just = "left",
                  gp = grid::gpar(fontfamily = "mono", cex = 0.65),
                  default.units = "npc")
  grid::popViewport()
  
  grid::popViewport()
}

model.info <- mod.info(final_model, Train.data, FALSE, FALSE)
cm_train <- ct_plot(
  model.info$class.table,
  plot.title = "Train Set (for model report)",
  conformation = "1. 1st Place"
)
grid_side_by_side_confusion(
  confusion_plot = cm_train$plot,
  confusion_table = model.info$class.table,
  left_title = "Train confusion matrix",
  right_title = "Train confusion table"
)

p_train_hm <- prob.heatmap(
  final_model,
  Train.data,
  plot.title = "Train Set (for model report)",
  conformation = "1. 1st Place"
)
print(p_train_hm)

if (nrow(Test.data) > 0) {
  model.info <- mod.info(final_model, Test.data, FALSE, FALSE)
  cm_test <- ct_plot(
    model.info$class.table,
    plot.title = "Test Set (for model report)",
    conformation = "1. 1st Place"
  )
  
  grid_side_by_side_confusion(
    confusion_plot = cm_test$plot,
    confusion_table = model.info$class.table,
    left_title = "Test confusion matrix",
    right_title = "Test confusion table"
  )

  p_test_hm <- prob.heatmap(
    final_model,
    Test.data,
    plot.title = "Test Set (for model report)",
    conformation = "1. 1st Place"
  )
  print(p_test_hm)
}

dev.off()
log_msg("Cross-validation model report PDF written successfully.")

log_msg("Finished Cross validation_Deuteration.R.")