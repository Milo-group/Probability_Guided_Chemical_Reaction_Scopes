library('rxn.cond.class')
library(caret)

# Train and test function to add to the McFadden values -------------------

fit_models_with_accuracy <- function(models_df, train_data, test_data) {
  
  n <- nrow(models_df)
  
  # Preallocate vectors
  train_acc <- numeric(n)
  test_acc  <- numeric(n)
  
  for (i in 1:n) {
    # Extract formula from first column
    test.form <- models_df[i, 1][[1]]
    
    # Fit multinomial model
    "For non-ordinal"
    #    model <- nnet::multinom(test.form, data = train_data, maxit = 2000, trace = FALSE) 
    
    "For ordinal"
    # Define starting coefficients
    num.of.vars <- stringi::stri_count(test.form, fixed = '+')
    start <- c(rep(0, num.of.vars + 2), 1)
    
    # Train model
    model <- fit_polr(formula = test.form, data = train_data)
    
    # Capture output from mod.info silently
    model.info.train <- { tmp <- capture.output(res <- mod.info(model, train_data, TRUE, TRUE)); res }
    model.info.test  <- { tmp <- capture.output(res <- mod.info(model, test_data, FALSE, FALSE)); res }
    
    # Store accuracies
    train_acc[i] <- model.info.train$accuracy_print
    test_acc[i]  <- round(model.info.test$accuracy, 2)
  }
  
  # Combine into final table
  final_table <- cbind(models_df, 'Train Accuracy' = train_acc, 'Test Accuracy' = test_acc)
  
  return(final_table)
}

# Train-test check --------------------------------------------------------

# Load data from a CSV file into a data frame
data <- data.frame(data.table::fread('Training_Data.csv'), check.names = F)

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

library(geometry)  # for convhulln
set.seed(1111)

# ------------#
# PCA
# ------------#

pca_res <- prcomp(Train.data[,2:ncol(Train.data)-1], center = TRUE, scale. = TRUE)

# Check importance
pca_var_percent <- summary(pca_res)$importance[2, ] * 100  # % variance per PC
cum_var_percent <- summary(pca_res)$importance[3, ] * 100  # cumulative %
round(cbind(pca_var_percent, cum_var_percent), 1)

"
     pca_var_percent cum_var_percent
PC1             28.4            28.4
PC2             22.9            51.3
PC3             14.6            65.9
PC4             11.8            77.6
PC5              7.3            85.0
PC6              4.0            89.0
PC7              3.2            92.2
...
"

core_train_pcs <- pca_res$x[, 1:3]  # first 3 PCs

candidate_test_pcs <- predict(pca_res, newdata = Test.data)[, 1:3]


# --------------#
# 3D convex hull
# --------------#

hull_indices <- convhulln(core_train_pcs, options = "FA")

inside_hull_test <- inhulln(hull_indices, candidate_test_pcs)

inhull_test_data <- Test.data[inside_hull_test, ] #in convex-hull
extreme_data <- Test.data[!inside_hull_test, ]

inhull_test_data <- Test.data[rownames(inhull_test_data), , drop = FALSE]
extreme_data <- Test.data[rownames(extreme_data), , drop = FALSE]

models <- sub_model_log(data = Train.data, 
                        min = 4, 
                        max = 4, 
                        ordinal = T)

knitr::kable(fit_models_with_accuracy(models, Train.data, inhull_test_data))

knitr::kable(fit_models_with_accuracy(models, Train.data, extreme_data))

"
Combination of the two tables gives the following table:

|formula                                                 | McFadden R2| Train Accuracy| Test Accuracy|  Extreme Accuracy|
|:-------------------------------------------------------|-----------:|--------------:|-------------:|-----------------:|
|`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`    |       0.577|          90.48|         88.89|            100.00|
|`class` ~ `-2-3-` + `dip_y` + `NPA_1` + `NPA_3`         |       0.577|          85.71|        100.00|            100.00|
|`class` ~ `dip_y` + `NPA_3` + `diff 1 2` + `L`          |       0.570|          78.57|         88.89|            100.00|
|`class` ~ `dip_y` + `Dist(2, 7)` + `NPA_3` + `diff 1 2` |       0.567|          85.71|        100.00|            100.00|
|`class` ~ `-2-3-` + `dip_y` + `NPA_1` + `NPA_2`         |       0.557|          80.95|         88.89|            100.00|
|`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `diff 1 2` |       0.555|          85.71|         88.89|            100.00|
|`class` ~ `-2-7-` + `dip_y` + `NPA_1` + `NPA_3`         |       0.555|          83.33|        100.00|            100.00|
|`class` ~ `dip_y` + `NPA_1` + `NPA_3` + `diff 2 7`      |       0.553|          83.33|        100.00|            100.00|
|`class` ~ `dip_y` + `NPA_3` + `NPA_7` + `diff 1 2`      |       0.551|          83.33|        100.00|            100.00|
|`class` ~ `-2-3-` + `Total` + `NPA_1` + `NPA_3`         |       0.549|          85.71|         77.78|            100.00|
"

# External validation and predictions -------------------------------------

# ----------------#
# --Convex-hull---#
# ----------------#

# Load data from a CSV file into a data frame
data <- data.frame(data.table::fread('Training_Data.csv'), check.names = F)

# Clean and organize data
"For training"
row.names(data) <- data[,2]  # Set the second column as row names
data$class <- as.factor(data$class)  # Convert the 'class' column to a factor
data <- data[,-c(1:2)]  # Remove the first and second columns (name and tag)

ExternalVal <- data.frame(data.table::fread('External_Validation_Data.csv'), check.names = F)
Predictions <- data.frame(data.table::fread('Predicting_New_Substrates_Data.csv'), check.names = F)
Predictions$class <- c("1","3","1","2")

"For test"
# Clean and organize data
candidate_test <- rbind(ExternalVal, Predictions)
row.names(candidate_test) <- candidate_test[,1]  # Set the second column as row names
candidate_test$class <- as.factor(candidate_test$class)  # Convert the 'class' column to a factor
candidate_test <- candidate_test[,-1]  # Remove the first and second columns (name and tag)

library(geometry)  # for convhulln
set.seed(1000)

# ------------#
# PCA
# ------------#
alldata <- rbind(data[,2:ncol(data)],candidate_test)
data_scaled <- scale(alldata[,!names(alldata) %in% c("class","flag")])

train.data <- cbind(data[,1],data_scaled[1:55,])
colnames(train.data)[1] <- "flag"
test.data <- data_scaled[56:dim(data_scaled)[1],]
  
pca_res <- prcomp(train.data[,2:ncol(train.data)], center = TRUE, scale. = TRUE)

# Check importance
pca_var_percent <- summary(pca_res)$importance[2, ] * 100  # % variance per PC
cum_var_percent <- summary(pca_res)$importance[3, ] * 100  # cumulative %
round(cbind(pca_var_percent, cum_var_percent), 1)

"
     pca_var_percent cum_var_percent
PC1             30.3            30.3
PC2             22.4            52.7
PC3             13.5            66.1
PC4             11.9            78.0
PC5              6.7            84.7
PC6              4.2            88.9
PC7              3.3            92.2
...
"

core_train_pcs <- pca_res$x[, 1:3]  # first 3 PCs

candidate_test_pcs <- predict(pca_res, newdata = test.data)[, 1:3]


# --------------#
# 3D convex hull
# --------------#


library(geometry)
hull_indices <- convhulln(core_train_pcs, options = "FA")

inside_hull_test <- inhulln(hull_indices, candidate_test_pcs)

inhull_test_data <- test.data[inside_hull_test, ] #in convex-hull
extreme_data <- test.data[!inside_hull_test, ]

train.data <- data[rownames(train.data), , drop = FALSE] #only one containing flag
inhull_test_data <- alldata[rownames(inhull_test_data), , drop = FALSE]
extreme_data <- alldata[rownames(extreme_data), , drop = FALSE]

models <- sub_model_log(data = train.data, 
                        min = 4, 
                        max = 4, 
                        ordinal = T)

knitr::kable(fit_models_with_accuracy(models, train.data, inhull_test_data))

knitr::kable(fit_models_with_accuracy(models, train.data, extreme_data))

"
|formula                                                 | McFadden R2| Train Accuracy| Test Accuracy|  Extreme Accuracy|
|:-------------------------------------------------------|-----------:|--------------:|-------------:|-----------------:|
|`class` ~ `dip_y` + `NPA_3` + `diff 1 2` + `loc.B5`     |       0.633|          87.27|        100.00|            100.00|
|`class` ~ `dip_y` + `NPA_1` + `NPA_3` + `loc.B5`        |       0.631|          85.45|        100.00|             85.71|
|`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`    |       0.626|          89.09|        100.00|            100.00|
|`class` ~ `-2-3-` + `dip_y` + `NPA_3` + `diff 1 2`      |       0.626|          85.45|        100.00|            100.00|
|`class` ~ `dip_y` + `Dist(2, 7)` + `NPA_3` + `diff 1 2` |       0.615|          90.91|        100.00|             85.71|
|`class` ~ `Total` + `NPA_3` + `diff 1 2` + `loc.B5`     |       0.612|          87.27|         66.67|            100.00|
|`class` ~ `dip_y` + `Dist(2, 7)` + `NPA_1` + `NPA_3`    |       0.611|          89.09|        100.00|             71.43|
|`class` ~ `Total` + `NPA_1` + `NPA_3` + `loc.B5`        |       0.611|          87.27|         66.67|            100.00|
|`class` ~ `Total` + `diff 1 2` + `diff 2 3` + `loc.B5`  |       0.610|          85.45|        100.00|            100.00|
|`class` ~ `dip_y` + `NPA_3` + `diff 1 2` + `L`          |       0.607|          81.82|        100.00|            100.00|
"
