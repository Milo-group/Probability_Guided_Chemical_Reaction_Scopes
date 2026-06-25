# Getting Started Example ג€” model construction
# Mirrors Aldehyde-Deuteration/Data/Model Construction_Deuteration.R on a tiny dummy dataset.

library("rxn.cond.class")

get_script_dir <- function() {
  path <- tryCatch(
    normalizePath(sys.frame(1)$ofile, winslash = "/"),
    error = function(e) NA_character_
  )
  if (!is.na(path)) return(dirname(path))
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  if (nzchar(file_arg)) return(dirname(normalizePath(file_arg, winslash = "/")))
  "."
}

script_dir <- get_script_dir()
training_csv <- file.path(script_dir, "Training_Data.csv")
external_csv <- file.path(script_dir, "External_Validation_Data.csv")
predict_csv <- file.path(script_dir, "Predicting_New_Substrates_Data.csv")

cat("=== Step 1: Load training data ===\n")
data <- data.frame(data.table::fread(training_csv), check.names = FALSE)
row.names(data) <- data[, 2]
data$class <- as.factor(data$class)
data <- data[, -c(1:2)]

cat(sprintf("Loaded %d molecules with %d descriptor columns.\n", nrow(data), ncol(data) - 2))

cat("\n=== Step 2: Similarity-based train / test split ===\n")
# On this small dataset, default sampling would select every molecule.
# We cap sample.size so ~2/3 go to training and ~1/3 remain for testing.
one <- simi.sampler(data, 1, sample.size = 6)
two <- simi.sampler(data, 2, sample.size = 6)
three <- simi.sampler(data, 3, sample.size = 6)
one_three <- simi.sampler(data, 1, 3, sample.size = 3)
two_three <- simi.sampler(data, 2, 3, sample.size = 3)
similarties <- unique(c(union(one, one_three), union(two, two_three), three))

Train.data <- data[similarties, ]
Test.data <- data[-similarties, ]
cat(sprintf("Training set: %d molecules | Test set: %d molecules\n", nrow(Train.data), nrow(Test.data)))

cat("\n=== Step 3: Model search (small example: 2 descriptors per model) ===\n")
models.ordinal <- sub_model_log(
  data = Train.data,
  min = 2,
  max = 2,
  ordinal = TRUE
)
cat("Top ranked models:\n")
print(knitr::kable(head(models.ordinal, 3)))

test.form <- models.ordinal[1, 1][[1]]
cat("\nBest model formula:", gsub("`", "", test.form), "\n")

cat("\n=== Step 4: Fit best model ===\n")
test <- fit_polr(formula = test.form, data = Train.data)

cat("\n--- Training set performance ---\n")
model.info <- mod.info(test, Train.data, TRUE, TRUE)
print(ct_plot(model.info$class.table, plot.title = "Training Set", conformation = "1. 1st Place")$plot)
print(prob.heatmap(test, Train.data, plot.title = "Training Set", conformation = "1. 1st Place"))

cat("\n--- Test set performance ---\n")
model.info <- mod.info(test, Test.data, TRUE, FALSE)
print(ct_plot(model.info$class.table, plot.title = "Test Set", conformation = "1. 1st Place")$plot)
print(prob.heatmap(test, Test.data, plot.title = "Test Set", conformation = "1. 1st Place"))

cat("\n=== Step 5: External validation ===\n")
External <- data.frame(data.table::fread(external_csv), check.names = FALSE)
RN <- External[, 1]
External <- External[, -1, drop = FALSE]
External$class <- as.factor(External$class)
row.names(External) <- RN
colnames(External) <- colnames(Train.data[, 2:ncol(Train.data)])

model.info <- mod.info(test, External, TRUE, FALSE)
print(ct_plot(model.info$class.table, plot.title = "External Validation", conformation = "1. 1st Place")$plot)
print(prob.heatmap(test, External, plot.title = "External Validation", conformation = "1. 1st Place"))

cat("\n=== Step 6: Predict new substrates ===\n")
Prediction.set <- data.frame(data.table::fread(predict_csv), check.names = FALSE)
RN <- Prediction.set[, 1]
Prediction.set <- Prediction.set[, -1, drop = FALSE]
row.names(Prediction.set) <- RN
colnames(Prediction.set) <- colnames(Train.data[, 2:(ncol(Train.data) - 1)])

pred_probs <- predict(test, Prediction.set, "probs") * 100
pred_class <- predict(test, Prediction.set, "class")
print(knitr::kable(cbind(pred_probs, predicted_class = pred_class)))

cat("\nDone. See Getting Started Example/README.md for what each step does.\n")
