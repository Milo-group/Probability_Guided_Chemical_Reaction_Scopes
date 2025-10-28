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

# Hydrogen Isotope Exchange -----------------------------------------------

# Import training data from CSV and convert it into a data frame
data50 <- data.frame(data.table::fread('HIE/Training_Data.csv')) 
row.names(data50) <- data50[,1]  # Set the first column as row names
data50 <- data50[,-1]  # Remove the first column after setting row names

# Convert all columns in the data frame to factors
for(i in 1:ncol(data50)) data50[, i] <- as.factor(data50[, i])

# Add a new column 'flag' that sequentially numbers the rows
data50 <- plyr::mutate(data50, flag = seq(1,nrow(data50)))

## Sampler data

# Import sampling data from CSV and convert to data frame
sampl.dat <- data.table::fread('HIE/Sampling_Data.csv') 
sampl.dat <- data.frame(sampl.dat)
row.names(sampl.dat) <- sampl.dat[,1]  # Set the first column as row names
sampl.dat <- sampl.dat[,-1]  # Remove the first column after setting row names

# Convert the 'class' column to a factor
sampl.dat$class <- as.factor(sampl.dat$class)

# Convert the first 60 columns of sampling data to numeric values
sampl.dat[,c(1:60)] <- as.numeric(as.matrix(sampl.dat[,c(1:60)]))

# Add a 'flag' column to sequentially number the rows
sampl.dat <- plyr::mutate(sampl.dat, flag = seq(1,nrow(sampl.dat)))

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

## Back to model search

# Combine similarities from the various groups sampled
similarties <- c(three,
                 union(two, two_three),
                 union(one, one_three),
                 union(four, four_three))

# Define training data from the subset of similarities
Train.data <- data50[similarties,]

# Define test data from the remaining rows not in the similarities set
Test.data <- data50[-similarties,]

test.form <- "`class` ~ `fr_NH0` + `fr_aldehyde` + `fr_amide` + `fr_pyridine`"

# 7-fold cross-validation
CV7 <- cv_classification(Train.data,test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)
round(CV7$overall_mean_accuracy*100,2)
"[1] 72.06"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)
round(LOO$overall_mean_accuracy*100,2)
"[1] 75"
