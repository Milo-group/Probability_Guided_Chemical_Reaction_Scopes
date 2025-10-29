library('rxn.cond.class')
library(caret)

# Deuteration -------------------------------------------------------------

# Load data from a CSV file into a data frame
data <- data.frame(data.table::fread('DeuterationData/Training_Data.csv'), check.names = F)

# Clean and organize data
row.names(data) <- data[,2]  # Set the second column as row names
data$class <- as.factor(data$class)  # Convert the 'class' column to a factor
data <- data[,-c(1:2)]  # Remove the first and second columns (name and tag)

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
# Extreme values ----------------------------------------------------------

# ----------- #
# --z-score-- #
# ----------- #

set.seed(1000)

# Compute z-scores for all descriptors
data_scaled <- scale(data[,!names(data) %in% c("class","flag")])

# Define extreme threshold (|z| > 3 in any descriptor)
extreme_idx <- which(apply(abs(data_scaled), 1, max) > 3)

" 2K  2R  2V 2AI 2AL 2AQ 2BA 2BB 2BC 
 11  18  22  35  38  43  53  54  55 "
# part data by extreme values to test (extreme) and train (normal).
extreme_data <- data[extreme_idx, ]
normal_data <- data[-extreme_idx, ]

models <- sub_model_log(data = normal_data, 
                                     min = 4, 
                                     max = 4, 
                                     ordinal = T)

fit_models_with_accuracy(models, normal_data, extreme_data)
"
|formula                                                 | McFadden R2| Train Accuracy| Test Accuracy|
|:-------------------------------------------------------|-----------:|--------------:|-------------:|
|`class` ~ `-2-3-` + `-2-7-` + `Total` + `Dist(1, 2)`    |       0.636|          86.96|         55.56|
|`class` ~ `-2-3-` + `-2-7-` + `dip_y` + `Dist(1, 2)`    |       0.624|          86.96|         55.56|
|`class` ~ `-2-3-` + `dip_y` + `Dist(1, 2)` + `diff 2 3` |       0.622|          84.78|         66.67|
|`class` ~ `dip_y` + `Dist(1, 2)` + `NPA_3` + `B5`       |       0.621|          82.61|         55.56|
|`class` ~ `dip_y` + `NPA_3` + `diff 1 2` + `B5`         |       0.613|          82.61|         77.78|
|`class` ~ `-2-3-` + `Total` + `Dist(1, 2)` + `diff 2 3` |       0.610|          86.96|         55.56|
|`class` ~ `dip_y` + `Total` + `Dist(1, 2)` + `NPA_3`    |       0.610|          82.61|         66.67|
|`class` ~ `-2-3-` + `dip_y` + `NPA_1` + `NPA_2`         |       0.609|          84.78|         77.78|
|`class` ~ `dip_y` + `Dist(1, 2)` + `NPA_3` + `L`        |       0.609|          84.78|         55.56|
|`class` ~ `-2-3-` + `dip_y` + `NPA_2` + `B5`            |       0.607|          84.78|         44.44|
"
# Preformance of chosen model in our original classification:

test.form <- "`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`"
num.of.vars <- stringi::stri_count(test.form, fixed = '+')
start <- c(rep(0, num.of.vars + 2), 1)
test <- fit_polr(formula = test.form, data = normal_data)

mod.info(test, normal_data, FALSE, FALSE)$McFadden
"0.59"
mod.info(test, normal_data, FALSE, FALSE)$accuracy
"84.78261"
mod.info(test, extreme_data, FALSE, FALSE)$accuracy
"77.77778"

# ------------ #
# -Electronic- #
# ------------ #
set.seed(1000)

# Compute z-scores for electronic descriptors
data_scaled <- scale(data[,names(data) %in% c("NPA_1","NPA_2","NPA_3","NPA_7")])

# Define extreme threshold (|z| > 3 in any descriptor)
extreme_idx <- which(apply(abs(data_scaled), 1, max) > 2)
"
2T  2V 2AZ 2BA 2BB 2BC 
 20  22  52  53  54  55
"
# part data by extreme values to test (extreme) and train (normal).
extreme_data <- data[extreme_idx, ]
normal_data <- data[-extreme_idx, ]

models <- sub_model_log(data = normal_data, 
                        min = 4, 
                        max = 4, 
                        ordinal = T)

knitr::kable(fit_models_with_accuracy(models, normal_data, extreme_data))
"
|formula                                                   | McFadden R2| Train Accuracy| Test Accuracy|
|:---------------------------------------------------------|-----------:|--------------:|-------------:|
|`class` ~ `-2-3-` + `-2-7-` + `Total` + `Dist(1, 2)`      |       0.615|          83.67|         33.33|
|`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`      |       0.567|          81.63|         66.67|
|`class` ~ `dip_y` + `NPA_3` + `diff 1 2` + `loc.B5`       |       0.562|          83.67|        100.00|
|`class` ~ `dip_x` + `dip_y` + `B5` + `loc.B5`             |       0.561|          85.71|         33.33|
|`class` ~ `-2-3-` + `dip_y` + `Dist(1, 2)` + `diff 2 3`   |       0.559|          81.63|         66.67|
|`class` ~ `-2-3-` + `dip_y` + `Dist(1, 2)` + `Dist(2, 7)` |       0.557|          81.63|         50.00|
|`class` ~ `-2-3-` + `dip_y` + `NPA_1` + `NPA_2`           |       0.557|          81.63|        100.00|
|`class` ~ `dip_y` + `NPA_1` + `NPA_3` + `loc.B5`          |       0.554|          75.51|        100.00|
|`class` ~ `-2-3-` + `dip_y` + `NPA_1` + `NPA_3`           |       0.552|          83.67|        100.00|
|`class` ~ `-2-3-` + `-2-7-` + `Dist(1, 2)` + `diff 2 3`   |       0.551|          75.51|         66.67|
"

# Preformance of chosen model in our original classification:

test.form <- "`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`"
num.of.vars <- stringi::stri_count(test.form, fixed = '+')
start <- c(rep(0, num.of.vars + 2), 1)
test <- fit_polr(formula = test.form, data = normal_data)

mod.info(test, normal_data, FALSE, FALSE)$McFadden
"0.567"
mod.info(test, normal_data, FALSE, FALSE)$accuracy
"81.63265"
mod.info(test, extreme_data, FALSE, FALSE)$accuracy
"66.66667"

# ---------- #
# --Steric-- #
# ---------- #

set.seed(1000)

# Compute z-scores for steric descriptors
data_scaled <- scale(data[,names(data) %in% c("B1","B5","L","loc.B5")])

# Define extreme threshold (|z| > 3 in any descriptor)
extreme_idx <- which(apply(abs(data_scaled), 1, max) > 2)
"
 2R  2U 2AI 2AL 2AO 2AQ 
 18  21  35  38  41  43 
"
# part data by extreme values to test (extreme) and train (normal).
extreme_data <- data[extreme_idx, ]
normal_data <- data[-extreme_idx, ]

models <- sub_model_log(data = normal_data, 
                        min = 4, 
                        max = 4, 
                        ordinal = T)

knitr::kable(fit_models_with_accuracy(models, normal_data, extreme_data))
"
|formula                                                 | McFadden R2| Train Accuracy| Test Accuracy|
|:-------------------------------------------------------|-----------:|--------------:|-------------:|
|`class` ~ `dip_y` + `NPA_3` + `diff 1 2` + `B5`         |       0.666|          89.80|         50.00|
|`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`    |       0.665|          89.80|         66.67|
|`class` ~ `dip_y` + `NPA_3` + `diff 1 2` + `loc.B5`     |       0.662|          87.76|         50.00|
|`class` ~ `dip_y` + `NPA_1` + `NPA_3` + `loc.B5`        |       0.652|          83.67|         50.00|
|`class` ~ `-2-3-` + `dip_y` + `NPA_2` + `diff 1 2`      |       0.645|          85.71|         66.67|
|`class` ~ `dip_y` + `Dist(2, 7)` + `NPA_1` + `NPA_3`    |       0.641|          85.71|         83.33|
|`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `diff 1 2` |       0.640|          83.67|         66.67|
|`class` ~ `dip_y` + `diff 1 2` + `diff 2 3` + `loc.B5`  |       0.640|          83.67|         50.00|
|`class` ~ `dip_y` + `NPA_1` + `NPA_3` + `diff 1 2`      |       0.636|          85.71|         83.33|
|`class` ~ `dip_y` + `NPA_1` + `NPA_3` + `B5`            |       0.634|          87.76|         66.67|"

# Preformance of chosen model in our original classification:

test.form <- "`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`"
num.of.vars <- stringi::stri_count(test.form, fixed = '+')
start <- c(rep(0, num.of.vars + 2), 1)
test <- fit_polr(formula = test.form, data = normal_data)

mod.info(test, normal_data, FALSE, FALSE)$McFadden
"0.665"
mod.info(test, normal_data, FALSE, FALSE)$accuracy
"89.79592"
mod.info(test, extreme_data, FALSE, FALSE)$accuracy
"66.66667"

# ----------------#
# --Convex-hull---#
# ----------------#

library(geometry)  # for convhulln
set.seed(1000)

# ------------#
# PCA
# ------------#
data_scaled <- scale(data[,!names(data) %in% c("class","flag")])

train_idx <- sample(1:nrow(data), size = round(0.75 * nrow(data)))
# PCA done on part of a random part of the data, since using the entire data would lead to no extremes
core_train <- data_scaled[train_idx, ]
candidate_test <- data_scaled[-train_idx, ]

pca_res <- prcomp(core_train, center = TRUE, scale. = TRUE)

# Check importance
pca_var_percent <- summary(pca_res)$importance[2, ] * 100  # % variance per PC
cum_var_percent <- summary(pca_res)$importance[3, ] * 100  # cumulative %
round(cbind(pca_var_percent, cum_var_percent), 1)

"
     pca_var_percent cum_var_percent
PC1             33.4            33.4
PC2             22.0            55.3
PC3             15.4            70.7
PC4              9.1            79.9
PC5              6.6            86.4
PC6              4.4            90.9
PC7              3.1            93.9
...
"

core_train_pcs <- pca_res$x[, 1:3]  # first 3 PCs

candidate_test_pcs <- predict(pca_res, newdata = candidate_test)[, 1:3]


# --------------#
# 3D convex hull
# --------------#


library(geometry)
hull_indices <- convhulln(core_train_pcs, options = "FA")

inside_hull_test <- inhulln(hull_indices, candidate_test_pcs)

train_data <- core_train
normal_test_data <- candidate_test[inside_hull_test, ] #in convex-hull
extreme_data <- candidate_test[!inside_hull_test, ]

train_data <- data[rownames(train_data), , drop = FALSE]
normal_test_data <- data[rownames(normal_test_data), , drop = FALSE]
extreme_data <- data[rownames(extreme_data), , drop = FALSE]

models <- sub_model_log(data = train_data, 
                        min = 4, 
                        max = 4, 
                        ordinal = T)

knitr::kable(fit_models_with_accuracy(models, train_data, normal_test_data))

knitr::kable(fit_models_with_accuracy(models, train_data, extreme_data))

"
|formula                                             | McFadden R2| Train Accuracy| Test Accuracy| Extreme Accuracy|
|:---------------------------------------------------|-----------:|--------------:|-------------:|----------------:|
|`class` ~ `dip_y` + `dip_z` + `diff 2 3` + `loc.B5` |       0.615|          82.93|         85.71|            71.43|
|`class` ~ `dip_y` + `NPA_3` + `diff 1 2` + `loc.B5` |       0.588|          82.93|        100.00|            71.43|
|`class` ~ `dip_y` + `dip_z` + `diff 2 7` + `loc.B5` |       0.586|          75.61|         71.43|            71.43|
|`class` ~ `dip_y` + `NPA_1` + `NPA_3` + `loc.B5`    |       0.571|          80.49|        100.00|            71.43|
|`class` ~ `-2-3-` + `dip_y` + `NPA_3` + `diff 1 2`  |       0.569|          85.37|        100.00|            71.43|
|`class` ~ `dip_y` + `dip_z` + `NPA_2` + `loc.B5`    |       0.562|          80.49|         85.71|            85.71|
|`class` ~ `-2-3-` + `dip_y` + `NPA_1` + `NPA_3`     |       0.559|          85.37|        100.00|            71.43|
|`class` ~ `-2-3-` + `dip_x` + `dip_y` + `NPA_1`     |       0.558|          80.49|        100.00|            85.71|
|`class` ~ `dip_y` + `NPA_3` + `NPA_7` + `diff 1 2`  |       0.558|          82.93|        100.00|            71.43|
|`class` ~ `dip_x` + `dip_y` + `diff 1 2` + `loc.B5` |       0.556|          80.49|         85.71|            85.71|
"

# Preformance of chosen model in our original classification:

test.form <- "`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`"
num.of.vars <- stringi::stri_count(test.form, fixed = '+')
start <- c(rep(0, num.of.vars + 2), 1)
test <- fit_polr(formula = test.form, data = train_data)

mod.info(test, train_data, FALSE, FALSE)$McFadden
"0.496"
mod.info(test, train_data, FALSE, FALSE)$accuracy
"80.4878"
mod.info(test, normal_test_data, FALSE, FALSE)$accuracy
"85.71429"
mod.info(test, extreme_data, FALSE, FALSE)$accuracy
"85.71429"

# --------------- #
# --SimiSampler-- #
# --------------- #

# Similarity-based sampling for each class
one <- simi.sampler(data, 1, plot=T)  # 2A, 2H extreme
two <- simi.sampler(data, 2, plot=T)  # 2AA, 2AF extreme
three <- simi.sampler(data, 3, plot=T)  # 2AU, 2BC extreme

# Define the training set from the subset of similarities
Train.data <- data[!row.names(data) %in% c("2A","2H","2AA","2AF","2AU","2BC"), ]

# Define the test set from the remaining rows not in the similarities set
Test.data <-  data[c("2A","2H","2AA","2AF","2AU","2BC"), ]

models.simi.sampler <- sub_model_log(data = Train.data, 
                         min = 4, 
                         max = 4, 
                         ordinal = T)

final_table1 <- fit_models_with_accuracy(models.simi.sampler, Train.data, Test.data)

knitr::kable(final_table1)
 "
|formula                                                | McFadden R2| Train Accuracy| Test Accuracy|
|:------------------------------------------------------|-----------:|--------------:|-------------:|
|`class` ~ `-2-3-` + `dip_y` + `NPA_3` + `diff 1 2`     |       0.779|          91.84|         50.00|
|`class` ~ `-2-3-` + `dip_y` + `NPA_1` + `diff 2 3`     |       0.763|          87.76|         66.67|
|`class` ~ `-2-3-` + `Total` + `NPA_1` + `NPA_3`        |       0.748|          89.80|         66.67|
|`class` ~ `dip_y` + `NPA_1` + `diff 1 2` + `diff 2 3`  |       0.744|          89.80|         50.00|
|`class` ~ `-2-3-` + `Total` + `NPA_3` + `diff 1 2`     |       0.742|          87.76|         66.67|
|`class` ~ `-2-3-` + `Total` + `diff 1 2` + `diff 2 3`  |       0.738|          85.71|         66.67|
|`class` ~ `Total` + `NPA_3` + `diff 1 2` + `diff 2 7`  |       0.735|          87.76|         66.67|
|`class` ~ `-2-3-` + `NPA_7` + `diff 1 2` + `diff 2 3`  |       0.733|          89.80|         66.67|
|`class` ~ `dip_y` + `NPA_3` + `diff 1 2` + `L`         |       0.733|          85.71|         50.00|
|`class` ~ `Total` + `diff 1 2` + `diff 2 3` + `loc.B5` |       0.728|          89.80|         50.00|
"
 # Preformance of chosen model in our original classification:
 
 test.form <- "`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`"
 num.of.vars <- stringi::stri_count(test.form, fixed = '+')
 start <- c(rep(0, num.of.vars + 2), 1)
 test <- fit_polr(formula = test.form, data = Train.data)
 
 mod.info(test, Train.data, FALSE, FALSE)$McFadden
 "0.662"
 mod.info(test, Train.data, FALSE, FALSE)$accuracy
 "85.71429"
 mod.info(test, Test.data, FALSE, FALSE)$accuracy
 "50"