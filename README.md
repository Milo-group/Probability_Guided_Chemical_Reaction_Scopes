# Probability Guided Chemical Reaction Scopes - Data Repository

This GitHub repository contains all data connected to our paper **"Probability Guided Chemical Reaction Scopes"**.  

The main folders are:

- **Aldehyde-Deuteration** ACS Catal. 2021, 11, 23, 14561–14569. <https://doi.org/10.1021/acscatal.1c04583>  
- **Hydrogen Isotope Exchange (HIE)**  ChemRxiv. 2023. <https://doi.org/10.26434/chemrxiv-2023-bmg14>  
- **Pyridine-Phosphination** J. Am. Chem. Soc. 2020, 142, 25, 11295–11305. <https://doi.org/10.1021/jacs.0c04674>  

These folders contain the data from the manuscript and all analyses that appear either in the manuscript or the Supporting Information (SI).  

Another subfolder is:

- **Extra case studies**  
  Contains data for other cases we modeled, which appear in the SI.  

---

## A guide for starters

If you are **not familiar with R or reaction-scope modeling**, begin with the folder **[Getting Started Example](Getting%20Started%20Example/)**. It contains:

- A **tiny dummy dataset** (made-up descriptor values, 30 molecules)
- Two scripts that mirror the real workflow in `Aldehyde-Deuteration/Data/`
- A short [README](Getting%20Started%20Example/README.md) with step-by-step instructions

**Quick start** (after installing `rxn.cond.class`):

1. In R, set your working directory to `Getting Started Example/Data`
2. Run `source("Model_Construction_Example.R")`, which loads data, finds a model, evaluates it, predicts new molecules
3. Optionally run `source("Cross_validation_Example.R")`, which checks whether accuracy is stable

When that makes sense, repeat the same steps on the real data in `Aldehyde-Deuteration/Data/`.

---

# General Guide for Using the Classification Code: Classify Chemical Reaction Conditions with the R Package `rxn.cond.class`

Below is a general guide for using the classification code provided in this repository.

`rxn.cond.class` is an R package designed to classify and visualize logistic regression classification models for chemical reaction conditions using both ordinal and non-ordinal models. It includes functionality for similarity-based sampling, model ranking, model evaluation, and heatmap visualization for model performance.

The workflow has four main stages: **install the package → prepare data and search models → fit and evaluate a model → predict on new substrates**. Each code block below corresponds to one step; comments marked **CHANGE** highlight settings you will typically adjust for your own dataset.

## Installation

### Installation from GitHub

Install `remotes` from CRAN first. The `repos = getCRANmirrors()[1,"URL"]` argument helps installation succeed in Linux interactive sessions.

```r
install.packages('remotes', repos = getCRANmirrors()[1,"URL"])
```

Install `rxn.cond.class` from GitHub, then load it. You only need to run the install line once per R installation.

```r
# Install (run once)
remotes::install_github('https://github.com/barkais/rxn.cond.class.git')

# Load (run at the start of every session)
library('rxn.cond.class')
```

Package overview and documentation are on the package [GitHub page](https://github.com/barkais/rxn.cond.class/tree/main).

For a hands-on walkthrough on a tiny dataset, see **[Getting Started Example](Getting%20Started%20Example/)** first.

## Example Usage

### Model Search

Load your training table, clean column names, split molecules into training and test sets using similarity-based sampling, and rank candidate models. The example below uses built-in package data; replace those lines with your own CSV when working on a real case study.

**Load and clean training data.** Row names should be molecule IDs; the `class` column must be a factor; descriptor columns remain as numeric predictors.

```r
# CHANGE: replace with your CSV path, e.g. data.table::fread("Training_Data.csv")
data <- rxn.cond.class::example_training_data

# Standard cleanup — adjust column indices if your CSV layout differs
row.names(data) <- data[, 2]          # CHANGE: column with molecule tag / ID
data$class <- as.factor(data$class)   # outcome must be a factor
data <- data[, -c(1:2)]               # CHANGE: drop name/tag columns (here: cols 1–2)
```

**Similarity-based train / test split.** For each class (or pair of classes), `simi.sampler()` picks representative molecules for training. On small datasets, pass `sample.size` so some molecules are left for testing (see Getting Started Example).

```r
# CHANGE: class labels (1, 2, 3) must match your data; add sample.size on small tables
one       <- simi.sampler(data, 1)           # class 1 vs itself
two       <- simi.sampler(data, 2)           # class 2 vs itself
three     <- simi.sampler(data, 3)           # class 3 vs itself
one_three <- simi.sampler(data, 1, 3)        # class 1 vs class 3
two_three <- simi.sampler(data, 2, 3)        # class 2 vs class 3

# Union of sampled row indices → training set
similarities <- c(union(one, one_three), union(two, two_three), three)
```

**Model search.** `sub_model_log()` ranks logistic models by accuracy. Set `min` / `max` to control how many descriptors each model may use; set `ordinal` to `TRUE` for ordered classes or `FALSE` for unordered classes.

```r
# CHANGE: min / max = number of descriptors per model; ordinal = TRUE or FALSE
models.ordinal <- sub_model_log(
  data = data[similarities, ],
  min = 2,           # CHANGE: minimum descriptors (e.g. 2)
  max = 2,           # CHANGE: maximum descriptors (e.g. 4 for a wider search)
  ordinal = TRUE
)

models.non.ordinal <- sub_model_log(
  data = data[similarities, ],
  min = 2,
  max = 2,
  ordinal = FALSE
)

# Training and held-out test sets from the similarity split
Train.set <- data[similarities, ]
Test.set  <- data[-similarities, ]
```

**Load external validation and prediction tables.** External data have known classes; prediction data do not. Match descriptor column names to the training set.

```r
# External validation (known labels) — CHANGE: your CSV or package example
External.set <- rxn.cond.class::example_validation_data
RN <- External.set$V1
External.set <- External.set[, -1]
External.set$class <- as.factor(External.set$class)
row.names(External.set) <- RN

# New substrates to predict (no labels) — CHANGE: your CSV or package example
Prediction.set <- rxn.cond.class::example_prediction_data
RN <- Prediction.set$V1
Prediction.set <- Prediction.set[, -1]
row.names(Prediction.set) <- RN
```

### Ordinal Model Example

Ordinal models treat class labels as ordered levels. Use this path when class order is chemically meaningful.

#### Model Ranking

Inspect the ranked list and pick a row. Row 1 is the top model by the package ranking criteria.

```r
knitr::kable(models.ordinal)
```

#### Training Set

Fit the top-ranked ordinal model with `fit_polr()` (proportional-odds logistic regression).

```r
# CHANGE: pick model rank (1 = best) or choose a formula from the table above
test.form <- models.ordinal[1, 1]

# Optional starting coefficients (usually leave as-is)
num.of.vars <- stringi::stri_count(test.form, fixed = '+')
start <- c(rep(0, num.of.vars + 2), 1)

# Fit on the training set
test <- fit_polr(formula = test.form, data = Train.set)
```

#### Model Information and Visualization (Training Set)

Summarize accuracy and plot the confusion matrix and probability heatmap on data the model was fit to.

```r
# mod.info: 3rd arg = print summary; 4th arg = include class probabilities
model.info <- mod.info(test, Train.set, TRUE, TRUE)

# CHANGE: plot.title and conformation are labels for your figures
confusion_matrix <- ct_plot(
  model.info$class.table,
  plot.title = 'Training Set',
  conformation = '1. 1st Place'   # CHANGE: short model label for the plot
)
confusion_matrix$plot

prob.heatmap(
  test, Train.set,
  plot.title = 'Training Set',
  conformation = '1. 1st Place'
)
```

#### Test Set

Evaluate on molecules excluded from training by the similarity split. This estimates performance on unseen substrates from the same dataset.

```r
# FALSE, FALSE = do not re-print full summary; evaluate on held-out rows only
model.info <- mod.info(test, Test.set, FALSE, FALSE)

confusion_matrix <- ct_plot(
  model.info$class.table,
  plot.title = 'Test Set',          # CHANGE: figure title
  conformation = '1. 1st Place'
)
confusion_matrix$plot

prob.heatmap(test, Test.set,
             plot.title = 'Test Set',
             conformation = '1. 1st Place')
```

#### External Validation

Evaluate on an independent dataset collected separately from the training table.

```r
model.info <- mod.info(test, External.set, FALSE)

confusion_matrix <- ct_plot(
  model.info$class.table,
  plot.title = 'External Validation',
  conformation = '1. 1st Place'
)
confusion_matrix$plot

prob.heatmap(test, External.set,
             plot.title = 'External Validation',
             conformation = '1. 1st Place')
```

#### Prediction of New Substrates

Predict class probabilities (%) and the most likely class for molecules without labels.

```r
knitr::kable(cbind(
  predict(test, Prediction.set, 'probs') * 100,   # class probabilities (%)
  predicted_class = predict(test, Prediction.set, 'class')
))
```

### Non-ordinal Model Example

Non-ordinal models treat classes as unordered categories. Use this path when class order is not meaningful or you want to compare against the ordinal approach.

#### Model Ranking

View the ranked non-ordinal models and choose which formula to fit.

```r
knitr::kable(models.non.ordinal)
```

#### Training Set

Fit a multinomial logistic regression with `nnet::multinom()`. Increase `maxit` if the optimizer fails to converge.

```r
# CHANGE: pick model rank (1 = best)
test.form <- models.non.ordinal[1, 1]

test <- nnet::multinom(
  test.form,
  data = Train.set,
  maxit = 2000,      # CHANGE: raise if you see convergence warnings
  trace = FALSE
)
```

#### Model Information and Visualization (Training Set)

Same plots as the ordinal workflow: confusion matrix and probability heatmap on the training set.

```r
model.info <- mod.info(test, Train.set, TRUE, TRUE)

confusion_matrix <- ct_plot(
  model.info$class.table,
  plot.title = 'Training Set',
  conformation = '1. 1st Place'
)
confusion_matrix$plot

prob.heatmap(test, Train.set,
             plot.title = 'Training Set',
             conformation = '1. 1st Place')
```

#### Test Set

Held-out performance for the non-ordinal model.

```r
model.info <- mod.info(test, Test.set, FALSE, FALSE)

confusion_matrix <- ct_plot(
  model.info$class.table,
  plot.title = 'Test Set',
  conformation = '1. 1st Place'
)
confusion_matrix$plot

prob.heatmap(test, Test.set,
             plot.title = 'Test Set',
             conformation = '1. 1st Place')
```

#### External Validation

Performance on the independent external validation set.

```r
model.info <- mod.info(test, External.set, FALSE)

confusion_matrix <- ct_plot(
  model.info$class.table,
  plot.title = 'External Validation',
  conformation = '1. 1st Place'
)
confusion_matrix$plot

prob.heatmap(test, External.set,
             plot.title = 'External Validation',
             conformation = '1. 1st Place')
```

## License

This package is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
