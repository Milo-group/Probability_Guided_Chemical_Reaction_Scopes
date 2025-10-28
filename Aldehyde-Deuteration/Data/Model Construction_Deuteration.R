"
# Install
remotes::install_github('https://github.com/barkais/rxn.cond.class.git')
"

library('rxn.cond.class')

# Deuteration -------------------------------------------------------------

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

# Model search with ordinal multinomial logistic regression
models.ordinal <- sub_model_log(data = Train.data, 
                                min = 4, 
                                max = 4, 
                                ordinal = T)

test.form <- models.ordinal[1,1]
# test.form <- "`class` ~ `-2-3-` + `dip_y` + `Dist(2, 7)` + `NPA_1`"

num.of.vars <- stringi::stri_count(test.form, fixed = '+')
start <- c(rep(0, num.of.vars + 2), 1)
test <- fit_polr(formula = test.form, data = Train.data)

# --------- #
# --Train-- #
# --------- #

# Confusion Matrix
model.info <- mod.info(test, Train.data, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'Training Set', 
        conformation = '1. 1st Place')$plot



# Heatmap
prob.heatmap(test, Train.data, 
             plot.title = 'Training Set', 
             conformation = '1. 1st Place')

# ---------- #
# -- Test -- #
# ---------- #

# Confusion Matrix
model.info <- mod.info(test, Test.data, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'Test Set', 
        conformation = '1. 1st Place')$plot

# Heatmap
prob.heatmap(test, Test.data, 
             plot.title = 'Test Set', 
             conformation = '1. 1st Place')

# ------------------------- #
# -- External Validation -- #
# ------------------------- #

External <- data.frame(data.table::fread('External_Validation_Data.csv'), check.names = F)
RN <- External$V1
External <- External[,-1]
External$class <- as.factor(External$class)
row.names(External) <- RN
colnames(External) <- colnames(Train.data[, c(2:dim(Train.data)[2])])

# Confusion Matrix
model.info <- mod.info(test, External, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'External Validation', 
        conformation = '1. 1st Place')$plot

# Heatmap
prob.heatmap(test, External, 
             plot.title = 'External Validation', 
             conformation = '1. 1st Place')

# ----------------- #
# -- Predictions -- #
# ----------------- #

# Load and organize prediction of new substrates data
Prediction.set <- data.frame(data.table::fread('Predicting_New_Substrates_Data.csv'), check.names = F)
RN <- Prediction.set$V1
Prediction.set <- Prediction.set[,-1]
row.names(Prediction.set) <- RN

knitr::kable(cbind(predict(test, Prediction.set, 'probs') * 100,
                   predicted_class = predict(test, Prediction.set, 'class')))
