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


# Halogenation through selective phosphination ----------------------------

#----------------#
#---  Binary ----#
#----------------#

# Import training data from a CSV file
data_bin <- data.frame(data.table::fread('Halogenation through phosphination/Binary/Training_Data_rdkit_stereoelectronic.csv'), check.names = F)

# Set the first column as row names
row.names(data_bin) <- data_bin[,1]  
data_bin <- data_bin[,-1]  # Remove the first column after setting row names

# Convert the 'class' column to a factor
data_bin$class <- as.factor(data_bin$class)

# Convert the first 50 columns to numeric values
data_bin[,c(1:50)] <- as.numeric(as.matrix(data_bin[,c(1:50)]))

# Convert columns 51 to 64 to factors after converting to numeric
for (i in 51:64) data_bin[,i] <- as.factor(as.numeric(data_bin[,i]))

# Random 10% out of group 1 ------------------------------------- #
set.seed(1000)

group1 <- which(data_bin$class == 1)

sample1 <- sample(group1, 5)  # 10% from group 1

Data_bin_noEx <- data_bin[-sample1, ]

ExternalValidation <- data_bin[sample1, ]

# Add a 'flag' column to sequentially number the rows
Data_bin_noEx <- plyr::mutate(Data_bin_noEx, flag = seq(1,nrow(Data_bin_noEx)))

# Sample from group 1 (first 50 columns + flag) with sample size being 75% of group 1 size
one <- simi.sampler(Data_bin_noEx[, c(1:50, 65:66)], 1, sample.size = round(sum(Data_bin_noEx$class == 1) * 0.75))

# Sample from group 2 (first 50 columns + flag) with sample size being 75% of group 2 size
two <- simi.sampler(Data_bin_noEx[, c(1:50, 65:66)], 2, sample.size = round(sum(Data_bin_noEx$class == 2) * 0.75))

# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
models <- sub_model_log(data = Data_bin_noEx[c(one, two), ], 
                        min = 4, 
                        max = 4, 
                        ordinal = F)

"|formula                                                                    | McFadden R2|
|:--------------------------------------------------------------------------|-----------:|
|`class` ~ `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene` + `fr_halogen`       |       0.721|
|`class` ~ `fr_aryl_methyl` + `fr_benzene` + `fr_halogen` + `fr_pyridine`   |       0.721|
|`class` ~ `NPA_1_P` + `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene`          |       0.689|
|`class` ~ `NPA_1_P` + `B1_2_SM` + `Iso.Polar_SM` + `fr_Ar_N`               |       0.683|
|`class` ~ `cross_P` + `fr_Ar_N` + `fr_aryl_methyl` + `fr_halogen`          |       0.676|
|`class` ~ `cross_P` + `fr_aryl_methyl` + `fr_halogen` + `fr_pyridine`      |       0.676|
|`class` ~ `NPA_1_P` + `Iso.Polar_SM` + `fr_Ar_N` + `fr_ether`              |       0.667|
|`class` ~ `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene` + `fr_ether`         |       0.662|
|`class` ~ `NPA_5_P` + `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene`          |       0.657|
|`class` ~ `NPA_Avg_H_onPhos` + `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene` |       0.644|"
 


# Define train and test data from the samples taken from groups 1 and 2
Train.data <- Data_bin_noEx[c(one, two), ]
Test.data <- Data_bin_noEx[-c(one, two), ]

# model 8 best in train-test

test.form <- "`class` ~ `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene` + `fr_ether`"

test <- nnet::multinom(test.form,
                       data = Data_bin_noEx[c(one, two), ],
                       maxit = 2000, 
                       trace = FALSE)

mod.info(test, Train.data, FALSE, FALSE)$accuracy
"93.02326"

mod.info(test, Test.data, FALSE, FALSE)$accuracy
"78.57143"
# External validation 

predict(test, newdata = ExternalValidation, type = "class")
"1 1 1 1 1
Levels: 1 2

Accuracy:100%"

CV7 <- cv_classification(Train.data[1:(ncol(Train.data)-1)], test.form, ordinal = FALSE,k = 7, 
                         n.iter = 200, seed = 100)

CV7$overall_mean_accuracy
"0.854898"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.data,test.form, ordinal = FALSE,k = "loo", 
                         n.iter = 1, seed = 100)
LOO$overall_mean_accuracy
"0.8837209"

# --------- #
# --Train-- #
# --------- #

# Confusion Matrix
model.info <- mod.info(test, Train.data, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'Training Set', 
        conformation = '8. 8th Place')$plot



# Heatmap
prob.heatmap(test, Train.data, 
             plot.title = 'Training Set', 
             conformation = '8. 8th Place')

# ---------- #
# -- Test -- #
# ---------- #

# Confusion Matrix
model.info <- mod.info(test, Test.data, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'Test Set', 
        conformation = '8. 8th Place')$plot

# Heatmap
prob.heatmap(test, Test.data, 
             plot.title = 'Test Set', 
             conformation = '8. 8th Place')

#---------------------#
#---  Multinomial ----#
#---------------------#

# Import training data from a CSV file
data_bin <- data.frame(data.table::fread("Halogenation through phosphination/3-class ordinal/Training_Data_rdkit_stereoelectronic.csv")) 

# Set the first column as row names
row.names(data_bin) <- data_bin[,1]  
data_bin <- data_bin[,-1]  # Remove the first column after setting row names

# Convert the 'class' column to a factor
data_bin$class <- as.factor(data_bin$class)  

# Convert the first 50 columns to numeric values
data_bin[,c(1:50)] <- as.numeric(as.matrix(data_bin[,c(1:50)]))  

# Convert columns 51 to 64 to factors after converting to numeric
for (i in 51:64) data_bin[,i] <- as.factor(as.numeric(data_bin[,i]))

# Random 10% out of group 2 ------------------------------------- #
set.seed(1000)

group2 <- which(data_bin$class == 2)

sample2 <- sample(group2, 5)  # 15% from group 2

Data_bin_noEx <- data_bin[-sample2, ]

# Add a 'flag' column to sequentially number the rows
Data_bin_noEx <- plyr::mutate(Data_bin_noEx, flag = seq(1,nrow(Data_bin_noEx)))

ExternalValidation <- data_bin[sample2, ]

# Sample from group 1 using simi.sampler
one <- simi.sampler(Data_bin_noEx, 1)

# Sample 21 molecules from group 2
two <- simi.sampler(Data_bin_noEx, 2, sample.size = 18)

# Sample from group 3 using simi.sampler
three <- simi.sampler(Data_bin_noEx, 3)

# Combine the similarities from groups 1, 2, and 3 into one vector
uni_similarties <- c(one, 
                     two, 
                     three)

models <- sub_model_log(data = Data_bin_noEx[uni_similarties, ], 
                        min = 3, 
                        max = 3, 
                        ordinal = T)

"                                                      formula McFadden R2
1          `class` ~ `para_P` + `NPA_6_SM` + `fr_aryl_methyl`       0.638
2         `class` ~ `NPA_1_P` + `loc.B5_1_SM` + `fr_bicyclic`       0.616
3           `class` ~ `Dist_.P.C._P` + `NPA_6_SM` + `para_SM`       0.607
4         `class` ~ `NPA_6_SM` + `para_SM` + `fr_aryl_methyl`       0.607
5     `class` ~ `NPA_5_P` + `Dist_.P.C._P` + `fr_aryl_methyl`       0.601
6             `class` ~ `NPA_1_P` + `Total_P` + `fr_bicyclic`       0.597
7  `class` ~ `Dist_.P.C._P` + `loc.B1_2_SM` + `para.angle_SM`       0.573
8              `class` ~ `para_P` + `NPA_6_SM` + `NPA_sum_SM`       0.569
9             `class` ~ `NPA_6_SM` + `NPA_sum_SM` + `para_SM`       0.568
10            `class` ~ `NPA_1_P` + `dip_y_P` + `fr_bicyclic`       0.563"

# Define training data from the subset of similarities
Train.Data <- Data_bin_noEx[uni_similarties, ]

# Define test data from the remaining rows not in the similarities set
Test.Data <- Data_bin_noEx[-uni_similarties, ]

# Train models on the subset of data corresponding to the combined similarities
test.form <- "`class` ~ `Dist_.P.C._P` + `NPA_6_SM` + `para_SM`"

# Define starting coefficients
num.of.vars <- stringi::stri_count(test.form, fixed = '+')
start <- c(rep(0, num.of.vars + 2), 1)

# Train model
test <- fit_polr(formula = test.form, data = Train.Data)

mod.info(test, Train.Data, FALSE, FALSE)$accuracy
"84.375"
mod.info(test, Test.Data, FALSE, FALSE)$accuracy
"100"

predict(test, newdata = ExternalValidation, type = "class")

"[1] 2 2 2 2 1
Levels: 1 2 3
Accuracy 80%"

CV7 <- cv_classification(Train.Data, test.form, ordinal = TRUE,k = 7, 
                         n.iter = 200, seed = 100)

CV7$overall_mean_accuracy
"0.8161071"

# Leave-one-out cross-validation
LOO <- cv_classification(Train.Data,test.form, ordinal = TRUE,k = "loo", 
                         n.iter = 1, seed = 100)
LOO$overall_mean_accuracy
"0.84375"

# --------- #
# --Train-- #
# --------- #

# Confusion Matrix
model.info <- mod.info(test, Train.Data, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'Training Set', 
        conformation = '3. 3rd Place')$plot



# Heatmap
prob.heatmap(test, Train.Data, 
             plot.title = 'Training Set', 
             conformation = '3. 3rd Place')

# ---------- #
# -- Test -- #
# ---------- #

# Confusion Matrix
model.info <- mod.info(test, Test.Data, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'Test Set', 
        conformation = '3. 3rd Place')$plot

# Heatmap
prob.heatmap(test, Test.Data, 
             plot.title = 'Test Set', 
             conformation = '3. 3rd Place')
