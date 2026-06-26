# Generate Model Search PDF reports from pre-ranked model tables
# (no sub_model_log / exhaustive search ג€” uses results recorded in McNally_10_percent_random.R)

library("rxn.cond.class")

get_script_dir <- function() {
  path <- tryCatch(
    normalizePath(sys.frame(1)$ofile, winslash = "/"),
    error = function(e) NA_character_
  )
  if (!is.na(path)) return(dirname(path))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
  if (nzchar(file_arg)) return(dirname(normalizePath(file_arg, winslash = "/")))
  "."
}

script_dir <- get_script_dir()
report_dir <- file.path(script_dir, "Model Reports")

source(file.path(script_dir, "model_search_report_helpers.R"), local = TRUE)

# Ranked models from McNally_10_percent_random.R (binary section) ----------------
models_binary <- data.frame(
  formula = c(
    "`class` ~ `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene` + `fr_halogen`",
    "`class` ~ `fr_aryl_methyl` + `fr_benzene` + `fr_halogen` + `fr_pyridine`",
    "`class` ~ `NPA_1_P` + `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene`",
    "`class` ~ `NPA_1_P` + `B1_2_SM` + `Iso.Polar_SM` + `fr_Ar_N`",
    "`class` ~ `cross_P` + `fr_Ar_N` + `fr_aryl_methyl` + `fr_halogen`",
    "`class` ~ `cross_P` + `fr_aryl_methyl` + `fr_halogen` + `fr_pyridine`",
    "`class` ~ `NPA_1_P` + `Iso.Polar_SM` + `fr_Ar_N` + `fr_ether`",
    "`class` ~ `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene` + `fr_ether`",
    "`class` ~ `NPA_5_P` + `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene`",
    "`class` ~ `NPA_Avg_H_onPhos` + `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene`"
  ),
  `McFadden R2` = c(0.721, 0.721, 0.689, 0.683, 0.676, 0.676, 0.667, 0.662, 0.657, 0.644),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Ranked models from McNally_10_percent_random.R (3-class ordinal section) -----
models_ordinal <- data.frame(
  formula = c(
    "`class` ~ `para_P` + `NPA_6_SM` + `fr_aryl_methyl`",
    "`class` ~ `NPA_1_P` + `loc.B5_1_SM` + `fr_bicyclic`",
    "`class` ~ `Dist_.P.C._P` + `NPA_6_SM` + `para_SM`",
    "`class` ~ `NPA_6_SM` + `para_SM` + `fr_aryl_methyl`",
    "`class` ~ `NPA_5_P` + `Dist_.P.C._P` + `fr_aryl_methyl`",
    "`class` ~ `NPA_1_P` + `Total_P` + `fr_bicyclic`",
    "`class` ~ `Dist_.P.C._P` + `loc.B1_2_SM` + `para.angle_SM`",
    "`class` ~ `para_P` + `NPA_6_SM` + `NPA_sum_SM`",
    "`class` ~ `NPA_6_SM` + `NPA_sum_SM` + `para_SM`",
    "`class` ~ `NPA_1_P` + `dip_y_P` + `fr_bicyclic`"
  ),
  `McFadden R2` = c(0.638, 0.616, 0.607, 0.607, 0.601, 0.597, 0.573, 0.569, 0.568, 0.563),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

prepare_binary_holdout_data <- function() {
  binary_csv <- file.path(script_dir, "Data/Binary/Training_Data_rdkit_stereoelectronic.csv")
  data_bin <- data.frame(data.table::fread(binary_csv), check.names = FALSE)
  row.names(data_bin) <- data_bin[, 1]
  data_bin <- data_bin[, -1]
  data_bin$class <- as.factor(data_bin$class)
  data_bin[, c(1:64)] <- as.numeric(as.matrix(data_bin[, c(1:64)]))

  set.seed(1000)
  sample1 <- sample(which(data_bin$class == 1), 5)
  Data_bin_noEx <- data_bin[-sample1, ]
  Data_bin_noEx <- plyr::mutate(Data_bin_noEx, flag = seq(1, nrow(Data_bin_noEx)))

  one <- simi.sampler(
    Data_bin_noEx[, c(1:50, 65:66)], 1,
    sample.size = round(sum(Data_bin_noEx$class == 1) * 0.75)
  )
  two <- simi.sampler(
    Data_bin_noEx[, c(1:50, 65:66)], 2,
    sample.size = round(sum(Data_bin_noEx$class == 2) * 0.75)
  )

  list(
    train = Data_bin_noEx[c(one, two), ],
    test = Data_bin_noEx[-c(one, two), ]
  )
}

prepare_ordinal_holdout_data <- function() {
  multi_csv <- file.path(script_dir, "Data/3-class ordinal/Training_Data_rdkit_stereoelectronic.csv")
  data_bin <- data.frame(data.table::fread(multi_csv))
  row.names(data_bin) <- data_bin[, 1]
  data_bin <- data_bin[, -1]
  data_bin$class <- as.factor(data_bin$class)
  data_bin[, c(1:64)] <- as.numeric(as.matrix(data_bin[, c(1:64)]))

  set.seed(1000)
  sample2 <- sample(which(data_bin$class == 2), 5)
  Data_bin_noEx <- data_bin[-sample2, ]
  Data_bin_noEx <- plyr::mutate(Data_bin_noEx, flag = seq(1, nrow(Data_bin_noEx)))

  one <- simi.sampler(Data_bin_noEx, 1)
  two <- simi.sampler(Data_bin_noEx, 2, sample.size = 18)
  three <- simi.sampler(Data_bin_noEx, 3)
  uni_similarties <- c(one, two, three)

  list(
    train = Data_bin_noEx[uni_similarties, ],
    test = Data_bin_noEx[-uni_similarties, ]
  )
}

log_msg("=== Generating Model Search reports from pre-ranked tables ===")

log_msg("Binary: preparing train/test split ...")
bin_split <- prepare_binary_holdout_data()
log_msg(sprintf("Binary: n=%d train, n=%d test", nrow(bin_split$train), nrow(bin_split$test)))

write_model_search_report(
  train_data = bin_split$train,
  test_data = bin_split$test,
  models = models_binary,
  ordinal = FALSE,
  report_title = "Pyridine Phosphination \u2013 Binary - Model Testing (10% holdout)",
  output_file = file.path(report_dir, "Model.Search.Pyr.Phos.Binary.10pctHoldout.pdf")
)

log_msg("3-class ordinal: preparing train/test split ...")
ord_split <- prepare_ordinal_holdout_data()
log_msg(sprintf("3-class ordinal: n=%d train, n=%d test", nrow(ord_split$train), nrow(ord_split$test)))

write_model_search_report(
  train_data = ord_split$train,
  test_data = ord_split$test,
  models = models_ordinal,
  ordinal = TRUE,
  report_title = "Pyridine Phosphination \u2013 3 categories - Model Testing (10% holdout)",
  output_file = file.path(report_dir, "Model.Search.Pyr.Phos.Ordinal.10pctHoldout.pdf")
)

log_msg("=== Done. Reports written to Model Reports/ ===")
