# Getting Started Example

## What is in this folder?

| File | Purpose |
|------|---------|
| `Data/Training_Data.csv` | Small training table (30 molecules, 3 classes, 6 descriptors) |
| `Data/External_Validation_Data.csv` | 3 held-out molecules with known class labels |
| `Data/Predicting_New_Substrates_Data.csv` | 2 new molecules **without** labels (for prediction only) |
| `Data/Model_Construction_Example.R` | Full workflow: split → search → fit → evaluate → predict |
| `Data/Cross_validation_Example.R` | Same split, then cross-validation checks |

The descriptor columns (`feat1` … `feat6`) do not represent real molecular properties. They contain dummy values designed to provide a simple demonstration in which the workflow can easily reach 100% performance.

On this small table, the scripts pass `sample.size` to `simi.sampler()` so some molecules are left out for testing. The full deuteration dataset in `Aldehyde-Deuteration/Data/` is large enough that this is not needed.

## Before you start

1. Install R (≥ 4.0 recommended).
2. Install the package (once):

```r
install.packages("remotes")
remotes::install_github("barkais/rxn.cond.class")
```

3. Install helper packages used by the scripts:

```r
install.packages(c("data.table", "knitr", "caret", "nnet"))
```

## How to run (3 steps)

### Step 1: Open R and set the working directory

In RStudio: **Session → Set Working Directory → Choose Directory…** and select the `Data` folder inside this example.

Or in the R console (edit the path to match your computer):

```r
setwd("path/to/Probability_Guided_Chemical_Reaction_Scopes/Getting Started Example/Data")
```

### Step 2: Build and evaluate a model

```r
source("Model_Construction_Example.R")
```

This script will:

1. Load `Training_Data.csv`
2. Split molecules into training and test sets using similarity sampling
3. Search for a good 2-descriptor model
4. Show accuracy tables and plots for train / test / external validation
5. Predict classes for the new substrates in `Predicting_New_Substrates_Data.csv`

### Step 3: Check stability with cross-validation (optional)

```r
source("Cross_validation_Example.R")
```

This runs 5-fold cross-validation (10 repeats) and leave-one-out validation on the training subset.

## What each CSV column means

| Column | Meaning |
|--------|---------|
| First column (empty header) | Molecule name / ID |
| `tag` | Short label (used as row name in R) |
| `flag` | Row index used by similarity sampling |
| `feat1` … `feat6` | Numeric descriptors (stand in for computed molecular features) |
| `class` | Outcome category (1, 2, or 3) |

## Moving to real data

Once you are comfortable with this example, replace the input CSV files with your own data and rerun the following scripts:

- `Model Construction_Deuteration.R`
- `Cross validation_Deuteration.R`

As long as the file names and table structure remain unchanged, the scripts should run without any additional modifications.
You can adjust the `sample.size` parameter (rows 33–38) to match the size of your dataset. You may specify either the actual number of samples or a percentage of the total dataset.

## Typical runtime

On the dummy data, each script should finish in under one minute.
