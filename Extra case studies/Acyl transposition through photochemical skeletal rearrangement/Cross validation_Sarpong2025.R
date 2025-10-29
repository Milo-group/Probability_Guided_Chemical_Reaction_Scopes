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

# Skeletal Editing --------------------------------------------------------

# Load data
Sarpong2025 <- read.csv('Data/computed_descriptors_Sarpong2025.csv')
rownames(Sarpong2025) <- Sarpong2025[,2]
Sarpong2025 <- Sarpong2025[,-c(1,2,3)]

# Convert the first 124 columns to numeric values
Sarpong2025[,c(1:124)] <- as.numeric(as.matrix(Sarpong2025[,c(1:124)]))

# data set as factor
for (i in 125:dim(Sarpong2025)[2]) Sarpong2025[,i] <- as.factor(as.numeric(Sarpong2025[,i]))

# Add a 'flag' column to sequentially number the rows
Sarpong2025 <- plyr::mutate(Sarpong2025, flag = seq(1,nrow(Sarpong2025)))

# Perform similarity-based sampling
# Sample from group 1 (first 124 columns + flag) with sample size being 75% of group 1 size
one <- simi.sampler(Sarpong2025[, c(1:124,(dim(Sarpong2025)[2]-1) ,dim(Sarpong2025)[2])], 1, sample.size = round(sum(Sarpong2025$class == 1) * 0.75))

# Sample from group 2 (first 124 columns + flag) with sample size being 75% of group 2 size
two <- simi.sampler(Sarpong2025[, c(1:124,(dim(Sarpong2025)[2]-1) ,dim(Sarpong2025)[2])], 2, sample.size = round(sum(Sarpong2025$class == 2) * 0.75))

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- Sarpong2025[c(one, two), 125:dim(Sarpong2025)[2]]
Test.data <- Sarpong2025[-c(one, two), 125:dim(Sarpong2025)[2]]

test.form <- "`class` ~ `fr_benzene` + `fr_methoxy` + `fr_para_hydroxylation`"

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

round(CV7$overall_mean_accuracy*100,2)
"[1] 78.23"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)
round(LOO$overall_mean_accuracy*100,2)
"[1] 77.42"
