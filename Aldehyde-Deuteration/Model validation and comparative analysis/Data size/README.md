# Progressive data-size evaluation

Measure how model performance changes as the training set grows.

## End-to-end workflow

### 1. Configure

Open `Progressive_DataSize_Evaluation.R` and edit the block at the bottom:

| Setting | Meaning |
|---------|---------|
| `split_method` | Train/holdout split: `random`, `similarity`, or `kennard_stone` |
| `initial_sampling` | How each subset is drawn: `random` or `stratified_min` |
| `min_per_class` | Minimum rows per class when `initial_sampling = "stratified_min"` |
| `balance_method` | Training-set balancing: `none`, `smote`, or `class_weight` |
| `start_size`, `max_size`, `step_add` | Subset sizes to evaluate (from first, to last, by step) |
| `train_p` | Holdout fraction when `split_method = "random"` |
| `n_repeats` | Independent repeats at each size |
| `n_candidates` | Formulas searched via `sub_model_log` on repeat 1 |
| `top_k_per_repeat` | Top formulas kept per repeat for stability ranking |
| `ordinal` | `TRUE` = ordered logit (`polr`), `FALSE` = multinomial |
| `seed` | Base random seed |

### 2. Run the analysis

From this folder:

```bash
Rscript Progressive_DataSize_Evaluation.R
```

Or in R:

```r
source("Progressive_DataSize_Evaluation.R")
```

Results are written to this folder unless you set `DATA_SIZE_OUTPUT_DIR`.

### 3. Visualize

Point the plotting script at the summary tables from step 2. One summary = one panel; several summaries = a comparison.

```bash
python data_size_accuracy_matplotlib.py \
  --summary <summary_1> <summary_2> \
  --labels "Condition A" "Condition B"
```

| Flag | Meaning |
|------|---------|
| `--summary` | Summary tables from the R run (one per condition) |
| `--raw` | Matching raw tables (inferred when omitted) |
| `--labels` | Panel / legend titles |
| `--data-dir` | Base folder for relative paths |
| `--output-dir` | Where figures are saved |
| `--out-prefix` | Prefix for output figure names |
| `--metrics` | Metrics for overlay figures when comparing 2+ series |

## Setup

**R:** `rxn.cond.class`, `caret`, `data.table`, `ggplot2`, `MASS`, `nnet`, `stringi`  
(add `recipes` and `themis` for SMOTE)

**Python:** `matplotlib`, `pandas`

**Data:** place training data under `<project>/Data/` or set `TRAINING_CSV`.  
Optional external validation: set `EXTERNAL_VALIDATION_CSV`.

## Files

| File | Purpose |
|------|---------|
| `Progressive_DataSize_Evaluation.R` | Run progressive data-size analysis |
| `data_size_accuracy_matplotlib.py` | Plot summary tables |
