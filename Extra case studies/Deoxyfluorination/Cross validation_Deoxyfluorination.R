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


# Deoxyfluorination 3-Cl sulfonyl ----------------------------------------

# Load data
ClSulfonyl <- read.csv('Data_noBase/Sulfonyl_3-Cl_Final.csv')
rownames(ClSulfonyl) <- ClSulfonyl[,2]
ClSulfonyl <- ClSulfonyl[,-c(1,2,82)] # Removes names and inchies

# Convert the first 158 columns to numeric values
ClSulfonyl[,c(17:79)] <- as.numeric(as.matrix(ClSulfonyl[,c(17:79)]))

# data set as factor
ClSulfonyl[, c(1:16, 80)] <- lapply(ClSulfonyl[, c(1:16, 80)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
ClSulfonyl <- plyr::mutate(ClSulfonyl, flag = seq(1,nrow(ClSulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(ClSulfonyl, 1)  # Sample from class 1
two <- simi.sampler(ClSulfonyl, 2)  # Sample from class 2
three <- simi.sampler(ClSulfonyl, 3) # Sample from class 3 
four <- simi.sampler(ClSulfonyl,4) # Sample from class 4

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(ClSulfonyl, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(ClSulfonyl, 2, 3) 

# Sample class 4 molecules similar to class 3
four_three <- simi.sampler(ClSulfonyl, 4, 3) 

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_three), 
                 union(two, two_three),
                 three,
                 union(four,four_three))

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- ClSulfonyl[similarties,]
Test.data <- ClSulfonyl[-similarties,]

test.form <- "`class` ~ `alcohol_.C1_exposed_area` + `ES_root_dipole` + `OC_L`"

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

round(CV7$overall_mean_accuracy*100,2)
"[1] 43.86"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)

round(LOO$overall_mean_accuracy*100,2)
"[1] 42.86"

# Deoxyfluorination Pyfluoro sulfonyl -----------------------------------------

# Load data
PyfluoroSulfonyl <- read.csv('Data_noBase/Sulfonyl_pyfluoro_Final.csv')
rownames(PyfluoroSulfonyl) <- PyfluoroSulfonyl[,2]
PyfluoroSulfonyl <- PyfluoroSulfonyl[,-c(1,2,82)] # Removes names and inchies

# Convert the first 158 columns to numeric values
PyfluoroSulfonyl[,c(17:79)] <- as.numeric(as.matrix(PyfluoroSulfonyl[,c(17:79)]))

# data set as factor
PyfluoroSulfonyl[, c(1:16, 80)] <- lapply(PyfluoroSulfonyl[, c(1:16, 80)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
PyfluoroSulfonyl <- plyr::mutate(PyfluoroSulfonyl, flag = seq(1,nrow(PyfluoroSulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(PyfluoroSulfonyl, 1)  # Sample from class 1
two <- simi.sampler(PyfluoroSulfonyl, 2,sample.size = round(0.7*sum(PyfluoroSulfonyl$class ==2)))  # Sample from class 2
three <- simi.sampler(PyfluoroSulfonyl, 3)  # Sample from class 3
four <- simi.sampler(PyfluoroSulfonyl, 4)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(PyfluoroSulfonyl, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(PyfluoroSulfonyl, 2, 3)

# Sample class 2 molecules similar to class 3
four_three <- simi.sampler(PyfluoroSulfonyl, 4, 3) 

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_three), 
                 union(two, two_three), 
                 three,
                 union(four,four_three))


# Define train and test data from the samples taken from groups 1 and 2
Train.data <- PyfluoroSulfonyl[similarties,]
Test.data <- PyfluoroSulfonyl[-similarties,]

test.form <- "`class` ~ `alcohol_.C1_exposed_area` + `alcohol_electronegativity` + `C_PVBur`"

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

round(CV7$overall_mean_accuracy*100,2)
"[1] 52.59"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)
round(LOO$overall_mean_accuracy*100,2)
"[1] 52.17"

# ============ #
# -- Binary -- #
# ============ #

# Load data
PyfluoroSulfonyl <- read.csv('Data_noBase/Sulfonyl_pyfluoro_Binary_NoBase_Final.csv')
rownames(PyfluoroSulfonyl) <- PyfluoroSulfonyl[,2]
PyfluoroSulfonyl <- PyfluoroSulfonyl[,-c(1,2,82)] # Removes names and inchies

# Convert the first 158 columns to numeric values
PyfluoroSulfonyl[,c(17:79)] <- as.numeric(as.matrix(PyfluoroSulfonyl[,c(17:79)]))

# data set as factor
PyfluoroSulfonyl[, c(1:16, 80)] <- lapply(PyfluoroSulfonyl[, c(1:16, 80)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
PyfluoroSulfonyl <- plyr::mutate(PyfluoroSulfonyl, flag = seq(1,nrow(PyfluoroSulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(PyfluoroSulfonyl, 1)  # Sample from class 1
two <- simi.sampler(PyfluoroSulfonyl, 2, sample.size = round(0.75*sum(PyfluoroSulfonyl$class == "2")))  # Sample from class 2

similarties <- c(one, two)

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- PyfluoroSulfonyl[similarties,]
Test.data <- PyfluoroSulfonyl[-similarties,]

# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
test.form <- "`class` ~ `alcohol_.C1_exposed_area` + `C_angle`"

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

round(CV7$overall_mean_accuracy*100,2)
"[1] 94.58"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)
round(LOO$overall_mean_accuracy*100,2)
"[1] 93.75"

# ================= #
# -- Multinomial -- #
# ================= #

# Load data
PyfluoroSulfonyl <- read.csv('Data_noBase/Sulfonyl_pyfluoro_Multi_NoBase_Final.csv')
rownames(PyfluoroSulfonyl) <- PyfluoroSulfonyl[,2]
PyfluoroSulfonyl <- PyfluoroSulfonyl[,-c(1,2,78)] # Removes names and inchies

# Convert the columns to numeric values
PyfluoroSulfonyl[,c(13:75)] <- as.numeric(as.matrix(PyfluoroSulfonyl[,c(13:75)]))

# data set as factor
PyfluoroSulfonyl[, c(1:12, 76)] <- lapply(PyfluoroSulfonyl[, c(1:12, 76)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
PyfluoroSulfonyl <- plyr::mutate(PyfluoroSulfonyl, flag = seq(1,nrow(PyfluoroSulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(PyfluoroSulfonyl, 1)  # Sample from class 1
two <- simi.sampler(PyfluoroSulfonyl, 2, sample.size = 7)  # Sample from class 2
three <- simi.sampler(PyfluoroSulfonyl, 3)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(PyfluoroSulfonyl, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(PyfluoroSulfonyl, 2, 3)

similarties <- c(union(one, one_three), union(two, two_three), three)

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- PyfluoroSulfonyl[similarties,]
Test.data <- PyfluoroSulfonyl[-similarties,]

# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
test.form <- "`class` ~ `alcohol_benzylic` + `C_PVBur`"

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

round(CV7$overall_mean_accuracy*100,2)
"[1] 52.41"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)

round(LOO$overall_mean_accuracy*100,2)
"[1] 44.44"

# Deoxyfluorination 3-CF3 sulfonyl ----------------------------------------

# Load data
CF3Sulfonyl <- read.csv('Data_noBase/Sulfonyl_3-CF3_Final.csv')
rownames(CF3Sulfonyl) <- CF3Sulfonyl[,2]
CF3Sulfonyl <- CF3Sulfonyl[,-c(1,2,82)] # Removes names and inchies

# Convert the first 158 columns to numeric values
CF3Sulfonyl[,c(17:79)] <- as.numeric(as.matrix(CF3Sulfonyl[,c(17:79)]))

# data set as factor
CF3Sulfonyl[, c(1:16, 80)] <- lapply(CF3Sulfonyl[, c(1:16, 80)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
CF3Sulfonyl <- plyr::mutate(CF3Sulfonyl, flag = seq(1,nrow(CF3Sulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(CF3Sulfonyl, 1)  # Sample from class 1
two <- simi.sampler(CF3Sulfonyl, 2)  # Sample from class 2
three <- simi.sampler(CF3Sulfonyl, 3)  # Sample from class 3
four <- simi.sampler(CF3Sulfonyl, 4)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(CF3Sulfonyl, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(CF3Sulfonyl, 2, 3)

# Sample class 2 molecules similar to class 3
four_three <- simi.sampler(CF3Sulfonyl, 4, 3) 

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_three), 
                 union(two, two_three), 
                 three,
                 union(four,four_three))


# Define train and test data from the samples taken from groups 1 and 2
Train.data <- CF3Sulfonyl[similarties,]
Test.data <- CF3Sulfonyl[-similarties,]
# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
test.form <- "`class` ~ `alcohol_.C1_exposed_area` + `ES_root_electronic_spatial_extent` + `ES_root_molar_volume`"

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

round(CV7$overall_mean_accuracy*100,2)
"[1] 54.78"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)

round(LOO$overall_mean_accuracy*100,2)
"[1] 52"

# Deoxyfluorination 3-NO2 sulfonyl -----------------------------------------

# Load data
NO2Sulfonyl <- read.csv('Data_noBase/Sulfonyl_3-NO2_Final.csv')
rownames(NO2Sulfonyl) <- NO2Sulfonyl[,2]
NO2Sulfonyl <- NO2Sulfonyl[,-c(1,2,82)] # Removes names and inchies

# Convert the first 158 columns to numeric values
NO2Sulfonyl[,c(17:79)] <- as.numeric(as.matrix(NO2Sulfonyl[,c(17:79)]))

# data set as factor
NO2Sulfonyl[, c(1:16, 80)] <- lapply(NO2Sulfonyl[, c(1:16, 80)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
NO2Sulfonyl <- plyr::mutate(NO2Sulfonyl, flag = seq(1,nrow(NO2Sulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(NO2Sulfonyl, 1,sample.size = round(0.45*sum(NO2Sulfonyl$class == 1)))  # Sample from class 1
two <- simi.sampler(NO2Sulfonyl, 2)  # Sample from class 2
three <- simi.sampler(NO2Sulfonyl, 3)  # Sample from class 3
four <- simi.sampler(NO2Sulfonyl, 4)  # Sample from class 3

# Sample class 1 molecules similar to class 2
one_two <- simi.sampler(NO2Sulfonyl, 1, 2)  

# Sample class 2 molecules similar to class 2
three_two <- simi.sampler(NO2Sulfonyl, 3, 2)

# Sample class 2 molecules similar to class 2
four_two <- simi.sampler(NO2Sulfonyl, 4, 2) 

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_two), 
                 two, 
                 union(three,three_two),
                 union(four,four_two))


# Define train and test data from the samples taken from groups 1 and 2
Train.data <- NO2Sulfonyl[similarties,]
Test.data <- NO2Sulfonyl[-similarties,]

test.form <- "`class` ~ `hardness` + `C_Mulliken_charge` + `O_ES_root_Mulliken_charge`"

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

round(CV7$overall_mean_accuracy*100,2)
"[1] 58.55"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)

round(LOO$overall_mean_accuracy*100,2)
"[1] 54.55"

# Deoxyfluorination PBSF sulfonyl -------------------------------------

# Load data
PBSFSulfonyl <- read.csv('Data_noBase/Sulfonyl_PBSF_Final.csv')
rownames(PBSFSulfonyl) <- PBSFSulfonyl[,2]
PBSFSulfonyl <- PBSFSulfonyl[,-c(1,2,81)] # Removes names and inchies

# Convert the first 158 columns to numeric values
PBSFSulfonyl[,c(16:78)] <- as.numeric(as.matrix(PBSFSulfonyl[,c(16:78)]))

# data set as factor
PBSFSulfonyl[, c(1:15, 79)] <- lapply(PBSFSulfonyl[, c(1:15, 79)], function(x) {
  as.factor(as.numeric(x))
})
# Add a 'flag' column to sequentially number the rows
PBSFSulfonyl <- plyr::mutate(PBSFSulfonyl, flag = seq(1,nrow(PBSFSulfonyl)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(PBSFSulfonyl, 1,sample.size = 0.5*sum(PBSFSulfonyl$class=="1"))  # Sample from class 1
two <- simi.sampler(PBSFSulfonyl, 2)  # Sample from class 2
three <- simi.sampler(PBSFSulfonyl, 3, sample.size = 4)  # Sample from class 3

# Sample class 1 molecules similar to class 2
one_two <- simi.sampler(PBSFSulfonyl, 1, 2)  

# Sample class 2 molecules similar to class 2
three_two <- simi.sampler(PBSFSulfonyl, 3, 2, sample.size = 4)

# Combine the similarities from the various classes into one vector

similarties <- c(union(one, one_two), 
                 two, 
                 union(three,three_two))

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- PBSFSulfonyl[similarties,]
Test.data <- PBSFSulfonyl[-similarties,]

test.form <- "`class` ~ `alcohol_secondary` + `ES_root_dipole` + `OC_L`"

CV7 <- cv_classification(Train.data, test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

round(CV7$overall_mean_accuracy*100,2)
"[1] 58.52"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)

round(LOO$overall_mean_accuracy*100,2)
"[1] 54.55"