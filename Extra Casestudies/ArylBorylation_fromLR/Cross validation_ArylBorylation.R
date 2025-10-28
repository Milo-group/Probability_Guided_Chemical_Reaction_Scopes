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

# Aryl Borylation ---------------------------------------------------------
# Load data
ArylBorylation <- read.csv('ArylBorylation_fromLR/Data/Doyle_data_organized_for_classification.csv')
rownames(ArylBorylation) <- ArylBorylation[,162]
ArylBorylation <- ArylBorylation[,-c(27,28,159:164)] # Removes names and smiles
ArylBorylation <- ArylBorylation[,-c(116,117)] #removel of all zero column

# Convert the first 158 columns to numeric values
ArylBorylation[,c(1:154)] <- as.numeric(as.matrix(ArylBorylation[,c(1:154)]))

# data set as factor
ArylBorylation[,155] <- as.factor(as.numeric(ArylBorylation[,155]))

# Add a 'flag' column to sequentially number the rows
ArylBorylation <- plyr::mutate(ArylBorylation, flag = seq(1,nrow(ArylBorylation)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(ArylBorylation, 1)  # Sample from class 1
two <- simi.sampler(ArylBorylation, 2)  # Sample from class 2
three <- simi.sampler(ArylBorylation, 3)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(ArylBorylation, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(ArylBorylation, 2, 3)  

# Combine the similarities from the various classes into one vector
"
similarties <- c(union(one, one_three), 
                 union(two, two_three), 
                 three)
                 "
similarties <- c(one,two,three)
# Define train and test data from the samples taken from groups 1 and 2
Train.data <- ArylBorylation[similarties,]
Test.data <- ArylBorylation[-similarties,]

test.form <- "`class` ~ `Kraken_E_solv_cds_boltz`"

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

round(CV7$overall_mean_accuracy*100,2)
"[1] 100"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)

round(LOO$overall_mean_accuracy*100,2)
"[1] 100"
