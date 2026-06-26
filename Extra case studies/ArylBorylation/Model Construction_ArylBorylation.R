"
# Install
remotes::install_github('barkais/rxn.cond.class', force = TRUE)
"

# Packege loading
library('rxn.cond.class')

# Load data
ArylBorylation <- read.csv(
  'Extra case studies/ArylBorylation/Data/Stevens_data_organized_for_classification_fixed.csv',
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(ArylBorylation) <- ArylBorylation[['Product Name']]

# Remove identifier / metadata columns by name
drop_cols <- c(
  'Electrophile', 'Electrophile_inchi',
  'Ligand', 'Ligand_inchi',
  'Product_inchi', 'Product Name old', 'Product Name',
  'Solvent_inchi', 'Yield'
)
ArylBorylation <- ArylBorylation[, !(names(ArylBorylation) %in% drop_cols)]

# Remove zero-variance numeric columns
is_zero_var <- vapply(
  ArylBorylation,
  function(x) is.numeric(x) && (length(unique(x)) == 1),
  logical(1)
)
ArylBorylation <- ArylBorylation[, !is_zero_var]

# Ensure numeric predictors and factor response
predictor_cols <- setdiff(names(ArylBorylation), 'class')
numeric_predictors <- lapply(
  ArylBorylation[predictor_cols],
  function(x) suppressWarnings(as.numeric(as.character(x)))
)
valid_predictors <- vapply(
  predictor_cols,
  function(col_name) {
    original <- ArylBorylation[[col_name]]
    converted <- numeric_predictors[[col_name]]
    all(is.na(converted) == is.na(original))
  },
  logical(1)
)

if (!all(valid_predictors)) {
  message('Dropping non-numeric predictors: ',
          paste(predictor_cols[!valid_predictors], collapse = ', '))
}

predictor_cols <- predictor_cols[valid_predictors]
ArylBorylation[predictor_cols] <- numeric_predictors[predictor_cols]
ArylBorylation$class <- suppressWarnings(as.numeric(as.character(ArylBorylation$class)))

# Keep only complete finite rows before similarity sampling
keep_rows <- complete.cases(ArylBorylation[, c(predictor_cols, 'class'), drop = FALSE]) &
  apply(ArylBorylation[, predictor_cols, drop = FALSE], 1, function(x) all(is.finite(x)))
ArylBorylation <- ArylBorylation[keep_rows, , drop = FALSE]

ArylBorylation$class <- as.factor(ArylBorylation$class)

# Add a 'flag' column to sequentially number the rows
ArylBorylation <- plyr::mutate(ArylBorylation, flag = seq(1,nrow(ArylBorylation)))

# Perform similarity-based sampling
# Similarity-based sampling for each class
one <- simi.sampler(ArylBorylation, 1, sample.size = 5)  # Sample from class 1
two <- simi.sampler(ArylBorylation, 2, sample.size = 5)  # Sample from class 2
three <- simi.sampler(ArylBorylation, 3, sample.size = 5)  # Sample from class 3

# Sample class 1 molecules similar to class 3
one_three <- simi.sampler(ArylBorylation, 1, 3)  

# Sample class 2 molecules similar to class 3
two_three <- simi.sampler(ArylBorylation, 2, 3)  

# Combine the similarities from the various classes into one vector
similarties <- c(one,two,three)

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- ArylBorylation[similarties,]
Test.data <- ArylBorylation[-similarties,]
# Train models using the McFadden approach on the subset of data corresponding to samples from groups 1 and 2
models <- sub_model_log(data = Train.data, 
                        min = 1, 
                        max = 3, 
                        ordinal = F)
knitr::kable(models)


# Training set ------------------------------------------------------------

# Use the first ranked non-ordinal model
test.form <- models[1, 1]

# Train the non-ordinal multinomial regression model
test <- nnet::multinom(test.form,
                       data = Train.data,
                       maxit = 2000, 
                       trace = FALSE)

# Cross-validation (smallest-group's-fold)
k.fold.log.iter(formula = test.form, 
                data = Train.data, 
                ordinal = FALSE, 
                stratify = TRUE, 
                iterations = 20, 
                verbose = TRUE)

# Leave-one-out cross-validation
k.fold.log.iter(formula = test.form, 
                data = Train.data, 
                ordinal = FALSE, 
                folds = nrow(Train.data), 
                stratify = FALSE, 
                iterations = 1, 
                verbose = TRUE)


# Visualization Training --------------------------------------------------

# Display model information and confusion matrix plot
model.info <- mod.info(test, Train.data, TRUE, TRUE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Training Set', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Train.data, 
             plot.title = 'Training Set', 
             conformation = '1. 1st Place')


# Test Set ----------------------------------------------------------------

# Evaluate the model on the test set
model.info <- mod.info(test, Test.data, FALSE, FALSE)

# Classification table plot
confusion_matrix <- ct_plot(model.info$class.table, 
                            plot.title = 'Test Set', 
                            conformation = '1. 1st Place')

confusion_matrix$plot

# Prediction probability heatmap
prob.heatmap(test, Test.data, 
             plot.title = 'Test Set', 
             conformation = '1. 1st Place')
