# Morgan feature tables for HIE modeling scripts (redundant bits removed).

.morgan_feature_spec <- function() {
  list(
    training_csv = "Training_morgan_filtered.csv",
    prediction_csv = "Prediction_morgan_filtered.csv",
    feature_pattern = "^MFP_[0-9]+$",
    label = "Morgan FP (redundant bits removed)",
    model_min = 4L,
    model_max = 4L,
    python_cmd = "python morgan_fp_filter.py"
  )
}

load_morgan_feature_table <- function(morgan_dir, spec = .morgan_feature_spec()) {
  training_csv <- file.path(morgan_dir, spec$training_csv)
  if (!file.exists(training_csv)) {
    stop(
      "Missing file:\n  ", training_csv,
      "\nRun ", spec$python_cmd, " in the Morgan FP folder."
    )
  }

  raw <- data.frame(data.table::fread(training_csv), check.names = FALSE)
  feat_cols <- grep(spec$feature_pattern, names(raw), value = TRUE)
  if (length(feat_cols) < 1L) {
    stop("No feature columns matched ", spec$feature_pattern, " in ", training_csv)
  }

  data50 <- raw
  row.names(data50) <- as.character(data50$Name)
  data50$Name <- NULL
  for (j in feat_cols) data50[[j]] <- as.numeric(data50[[j]])
  data50$class <- factor(data50$class)
  data50 <- plyr::mutate(data50, flag = seq_len(nrow(data50)))

  list(
    spec = spec,
    training_csv = training_csv,
    prediction_csv = file.path(morgan_dir, spec$prediction_csv),
    feature_cols = feat_cols,
    data = data50
  )
}
