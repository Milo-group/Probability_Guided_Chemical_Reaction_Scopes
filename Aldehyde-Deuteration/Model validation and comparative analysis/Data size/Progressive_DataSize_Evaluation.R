# Progressive data-size evaluation: grow training subsets, search formulas, export metrics.

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(file_path) > 0 && nzchar(file_path[1])) {
    return(dirname(normalizePath(file_path[1], winslash = "/", mustWork = FALSE)))
  }

  src_path <- NULL
  src_path <- tryCatch(
    normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE),
    error = function(e) NULL
  )
  if (!is.null(src_path) && nzchar(src_path)) {
    return(dirname(src_path))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()

ks_sampling_path <- normalizePath(
  file.path(script_dir, "..", "Kennard-stone", "kennard_stone_sampling.R"),
  winslash = "/",
  mustWork = FALSE
)

ensure_kennard_stone_helpers <- function(re_source = TRUE) {
  if (!file.exists(ks_sampling_path)) {
    stop(
      "Kennard–Stone helpers not found at: ", ks_sampling_path,
      "\nExpected: Model validation and comparative analysis/Kennard-stone/kennard_stone_sampling.R",
      call. = FALSE
    )
  }
  needs_load <- re_source ||
    !exists("sample_kennard_stone_indices", mode = "function") ||
    !exists("split_with_kennard_stone_sampler", mode = "function") ||
    !exists("prepare_training_data", mode = "function") &&
    !exists("prepare_deuteration_training_data", mode = "function")
  if (needs_load) {
    source(ks_sampling_path, local = FALSE)
    message("Loaded Kennard–Stone helpers: ", ks_sampling_path)
  }
  invisible(TRUE)
}

if (file.exists(ks_sampling_path)) {
  ensure_kennard_stone_helpers(re_source = TRUE)
} else {
  warning("Kennard-stone helpers not found at: ", ks_sampling_path)
}

find_package_dir <- function(start_dir, package_name = "rxn.cond.class", max_up = 10L) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 0:max_up) {
    candidate <- file.path(cur, package_name)
    if (dir.exists(candidate) && file.exists(file.path(candidate, "DESCRIPTION"))) {
      return(candidate)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  NULL
}

pkg_dir <- find_package_dir(script_dir)
if (!is.null(pkg_dir) && requireNamespace("devtools", quietly = TRUE)) {
  message(sprintf("Development mode: loading local package from %s", pkg_dir))
  devtools::load_all(pkg_dir, quiet = TRUE)
} else {
  library("rxn.cond.class")
}

split_with_random_partition <- function(data, seed = NULL, train_p = 0.75) {
  if (!requireNamespace("caret", quietly = TRUE)) stop("Please install 'caret'")
  if (!is.null(seed)) set.seed(seed)
  tr_idx <- caret::createDataPartition(data$class, p = train_p, list = FALSE)
  list(
    train_data = data[tr_idx, , drop = FALSE],
    test_data = data[-tr_idx, , drop = FALSE]
  )
}

split_with_similarity_sampler <- function(data, seed = 42, target_train_n = NULL) {
  if (!is.null(seed)) set.seed(seed)

  stratified_fallback <- function(df) {
    cls <- sort(unique(as.character(df$class)))
    picked <- integer(0)
    for (cl in cls) {
      idx <- which(as.character(df$class) == cl)
      if (length(idx) == 0) next
      take_n <- max(1, floor(length(idx) * 0.7))
      take_n <- min(take_n, length(idx))
      picked <- c(picked, sample(idx, size = take_n))
    }
    picked <- unique(picked)
    if (length(picked) == 0) {
      picked <- sample(seq_len(nrow(df)), size = max(1, floor(0.7 * nrow(df))))
    }
    list(
      train_data = df[picked, , drop = FALSE],
      test_data = df[-picked, , drop = FALSE]
    )
  }

  classes_present <- sort(unique(as.character(data$class)))
  if (!all(c("1", "2", "3") %in% classes_present)) {
    return(stratified_fallback(data))
  }

  similarities <- tryCatch(
    {
      one <- simi.sampler(data, 1)
      two <- simi.sampler(data, 2)
      three <- simi.sampler(data, 3)
      one_three <- simi.sampler(data, 1, 3)
      two_three <- simi.sampler(data, 2, 3)
      unique(c(union(one, one_three), union(two, two_three), three))
    },
    error = function(e) {
      message(sprintf(
        "simi.sampler failed on subset (n=%d): %s. Using stratified fallback split.",
        nrow(data), conditionMessage(e)
      ))
      NULL
    }
  )
  if (is.null(similarities) || length(similarities) == 0) {
    return(stratified_fallback(data))
  }

  similarities <- as.vector(similarities)
  similarities <- suppressWarnings(as.integer(similarities))
  similarities <- similarities[!is.na(similarities)]
  similarities <- unique(similarities)
  if ("flag" %in% names(data)) {
    idx_from_flag <- which(data$flag %in% similarities)
    idx_from_flag <- unique(idx_from_flag)
    if (length(idx_from_flag) == 0) {
      return(stratified_fallback(data))
    }
    similarities <- idx_from_flag
  }
  similarities <- similarities[similarities >= 1 & similarities <= nrow(data)]
  if (length(similarities) == 0) {
    return(stratified_fallback(data))
  }

  train_data <- data[similarities, , drop = FALSE]
  test_data <- data[-similarities, , drop = FALSE]
  if (nrow(train_data) == 0 || nrow(test_data) == 0) {
    return(stratified_fallback(data))
  }
  list(train_data = train_data, test_data = test_data)
}

align_class_levels <- function(df, class_levels) {
  if (!is.null(df) && "class" %in% names(df)) {
    df$class <- factor(as.character(df$class), levels = class_levels)
  }
  df
}

safe_accuracy <- function(model, eval_data, test.form) {
  if (is.null(eval_data) || nrow(eval_data) == 0) return(NA_real_)
  resp_name <- gsub("`", "", deparse(as.formula(test.form)[[2]]))
  pred <- tryCatch(
    predict(model, newdata = eval_data, type = "class"),
    error = function(e) rep(NA, nrow(eval_data))
  )
  pred <- as.vector(pred)
  truth <- as.vector(eval_data[[resp_name]])
  valid <- !is.na(pred) & !is.na(truth)
  if (!any(valid)) return(NA_real_)
  mean(as.character(pred[valid]) == as.character(truth[valid]))
}

safe_recall_for_class <- function(model, eval_data, test.form, target_class) {
  if (is.null(eval_data) || nrow(eval_data) == 0) return(NA_real_)
  resp_name <- gsub("`", "", deparse(as.formula(test.form)[[2]]))
  truth <- as.vector(eval_data[[resp_name]])
  if (is.null(truth)) return(NA_real_)
  truth_chr <- as.character(truth)
  target_idx <- which(truth_chr == as.character(target_class))
  if (length(target_idx) == 0) return(NA_real_)

  pred <- tryCatch(
    predict(model, newdata = eval_data, type = "class"),
    error = function(e) rep(NA, nrow(eval_data))
  )
  pred_chr <- as.character(as.vector(pred))
  mean(pred_chr[target_idx] == truth_chr[target_idx], na.rm = TRUE)
}

compute_evaluation_metrics <- function(model,
                                       train_data,
                                       holdout_test_data,
                                       all_data,
                                       sampled_idx,
                                       external_data,
                                       test.form,
                                       smallest_class) {
  if (is.null(model)) {
    return(list(
      train_accuracy = NA_real_,
      test_accuracy = NA_real_,
      external_accuracy = NA_real_,
      smallest_group_recall = NA_real_
    ))
  }
  external_eval_data <- build_external_eval_data(
    all_data = all_data,
    selected_idx = sampled_idx,
    external_data = external_data
  )
  list(
    train_accuracy = safe_accuracy(model, train_data, test.form),
    test_accuracy = safe_accuracy(model, holdout_test_data, test.form),
    external_accuracy = safe_accuracy(model, external_eval_data, test.form),
    smallest_group_recall = safe_recall_for_class(
      model = model,
      eval_data = all_data,
      test.form = test.form,
      target_class = smallest_class
    )
  )
}

mcfadden_col_name <- function(models_df) {
  nm <- names(models_df)[grep("McFadden", names(models_df), ignore.case = TRUE)]
  if (length(nm) == 0) NA_character_ else nm[1L]
}

compute_mcfadden_polr <- function(model, train_data, formula_chr, weights = NULL) {
  if (is.null(model)) return(NA_real_)
  null_form <- tryCatch(
    stats::as.formula(paste(deparse(stats::as.formula(formula_chr)[[2]]), "~ 1")),
    error = function(e) NULL
  )
  if (is.null(null_form)) return(NA_real_)
  m0_args <- list(formula = null_form, data = train_data, Hess = TRUE)
  if (!is.null(weights)) m0_args$weights <- weights
  m0 <- tryCatch(do.call(MASS::polr, m0_args), error = function(e) NULL)
  if (is.null(m0)) return(NA_real_)
  ll0 <- as.numeric(stats::logLik(m0))
  ll1 <- as.numeric(stats::logLik(model))
  if (!is.finite(ll0) || ll0 == 0 || !is.finite(ll1)) return(NA_real_)
  round(1 - (ll1 / ll0), digits = 3)
}

fit_model_from_formula <- function(formula_chr,
                                     train_data,
                                     ordinal = TRUE,
                                     weights = NULL) {
  formula_chr <- as.character(formula_chr)[1]
  if (!nzchar(formula_chr)) return(NULL)

  if (!ordinal) {
    form <- tryCatch(stats::as.formula(formula_chr), error = function(e) NULL)
    if (is.null(form)) return(NULL)
    args <- list(form, data = train_data, maxit = 2000, trace = FALSE)
    if (!is.null(weights)) args$weights <- weights
    return(tryCatch(do.call(nnet::multinom, args), error = function(e) NULL))
  }

  if (!is.null(weights)) {
    return(fit_polr_weighted(formula_chr, train_data, weights))
  }

  form <- tryCatch(stats::as.formula(formula_chr), error = function(e) NULL)
  if (is.null(form)) return(NULL)

  num.of.vars <- stringi::stri_count(formula_chr, fixed = "+")
  start <- c(rep(0, num.of.vars + 2), 1)
  success <- FALSE
  model <- NULL
  while (!success) {
    res <- tryCatch(
      list(
        ok = TRUE,
        model = MASS::polr(
          form,
          data = train_data,
          Hess = TRUE,
          start = start,
          control = list(maxit = 100)
        )
      ),
      error = function(e) list(ok = FALSE)
    )
    if (isTRUE(res$ok)) {
      success <- TRUE
      model <- res$model
    } else {
      start <- c(0, start)
      if (length(start) > num.of.vars + 15L) break
    }
  }

  if (is.null(model)) {
    model <- tryCatch(
      fit_polr(formula = formula_chr, data = train_data),
      error = function(e) NULL
    )
  }
  model
}

evaluate_formula_on_split <- function(formula_chr,
                                      train_data,
                                      holdout_test_data,
                                      all_data,
                                      sampled_idx,
                                      external_data,
                                      smallest_class,
                                      train_for_fit = NULL,
                                      case_weights = NULL,
                                      ordinal = TRUE) {
  if (is.null(train_for_fit)) train_for_fit <- train_data
  fit_weights <- case_weights
  if (!is.null(fit_weights)) {
    if (length(fit_weights) != nrow(train_for_fit)) {
      fit_weights <- make_class_weights_vec(train_for_fit$class)
    }
  }
  model <- fit_model_from_formula(
    formula_chr,
    train_for_fit,
    ordinal = ordinal,
    weights = fit_weights
  )
  mcf <- if (is.null(model) || !ordinal) {
    NA_real_
  } else {
    compute_mcfadden_polr(model, train_for_fit, formula_chr, weights = fit_weights)
  }
  met <- compute_evaluation_metrics(
    model = model,
    train_data = train_data,
    holdout_test_data = holdout_test_data,
    all_data = all_data,
    sampled_idx = sampled_idx,
    external_data = external_data,
    test.form = formula_chr,
    smallest_class = smallest_class
  )
  c(met, list(mcfadden_train = mcf))
}

evaluate_formula_candidates <- function(formulas,
                                        train_data,
                                        holdout_test_data,
                                        all_data,
                                        sampled_idx,
                                        external_data,
                                        smallest_class,
                                        train_for_fit = NULL,
                                        case_weights = NULL,
                                        mcfadden_train = NULL,
                                        ordinal = TRUE) {
  formulas <- unique(as.character(formulas))
  formulas <- formulas[nzchar(formulas)]
  if (length(formulas) == 0) {
    return(data.frame(
      formula = character(0),
      train_accuracy = numeric(0),
      test_accuracy = numeric(0),
      external_accuracy = numeric(0),
      smallest_group_recall = numeric(0),
      mcfadden_train = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  rows <- vector("list", length(formulas))
  for (i in seq_along(formulas)) {
    f <- formulas[[i]]
    met <- evaluate_formula_on_split(
      formula_chr = f,
      train_data = train_data,
      holdout_test_data = holdout_test_data,
      all_data = all_data,
      sampled_idx = sampled_idx,
      external_data = external_data,
      smallest_class = smallest_class,
      train_for_fit = train_for_fit,
      case_weights = case_weights,
      ordinal = ordinal
    )
    mcf <- met$mcfadden_train
    if (!is.null(mcfadden_train) && f %in% names(mcfadden_train) && is.na(mcf)) {
      mcf <- as.numeric(mcfadden_train[[f]])
    }
    rows[[i]] <- data.frame(
      formula = f,
      train_accuracy = met$train_accuracy,
      test_accuracy = met$test_accuracy,
      external_accuracy = met$external_accuracy,
      smallest_group_recall = met$smallest_group_recall,
      mcfadden_train = mcf,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

rank_top_formulas <- function(eval_df, top_k = 10L) {
  if (is.null(eval_df) || nrow(eval_df) == 0) return(character(0))
  top_k <- min(as.integer(top_k)[1L], nrow(eval_df))
  mcf <- eval_df$mcfadden_train
  acc <- eval_df$train_accuracy
  if (!all(is.na(mcf))) {
    ord <- order(mcf, acc, eval_df$test_accuracy, decreasing = TRUE, na.last = TRUE)
  } else if (!all(is.na(acc))) {
    ord <- order(acc, eval_df$test_accuracy, decreasing = TRUE, na.last = TRUE)
  } else {
    ord <- seq_len(nrow(eval_df))
  }
  ord <- ord[seq_len(min(top_k, length(ord)))]
  as.character(eval_df$formula[ord])
}

normalize_formula_chr <- function(x) {
  trimws(as.character(x))
}

summarize_top10_appearances <- function(top10_by_repeat) {
  all_formulas <- unique(unlist(top10_by_repeat, use.names = FALSE))
  if (length(all_formulas) == 0) {
    return(data.frame(
      formula = character(0),
      top10_count = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  counts <- vapply(all_formulas, function(f) {
    sum(vapply(top10_by_repeat, function(x) {
      if (is.null(x) || length(x) == 0) return(FALSE)
      f %in% x
    }, logical(1)))
  }, integer(1))
  out <- data.frame(
    formula = all_formulas,
    top10_count = as.integer(counts),
    stringsAsFactors = FALSE
  )
  out[order(-out$top10_count, out$formula), , drop = FALSE]
}

safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x)
}

BALANCE_METHODS <- c("none", "smote", "class_weight")
SPLIT_METHODS <- c("random", "similarity", "kennard_stone")
INITIAL_SAMPLING_METHODS <- c("random", "stratified_min")

balance_path_tag <- function(balance_method = BALANCE_METHODS) {
  balance_method <- match.arg(balance_method, BALANCE_METHODS)
  switch(
    balance_method,
    none = "",
    smote = "_smote",
    class_weight = "_classweight"
  )
}

progressive_csv_basename <- function(file_kind,
                                      split_method = SPLIT_METHODS,
                                      initial_sampling = INITIAL_SAMPLING_METHODS,
                                      balance_method = BALANCE_METHODS) {
  balance_method <- match.arg(balance_method, BALANCE_METHODS)
  split_method <- match.arg(split_method, SPLIT_METHODS)
  initial_sampling <- match.arg(initial_sampling, INITIAL_SAMPLING_METHODS)
  if (balance_method == "none") {
    tag <- switch(split_method, random = "random", similarity = "simi", kennard_stone = "ks")
    start_tag <- switch(
      initial_sampling,
      stratified_min = "stratstart",
      "randomstart"
    )
    return(sprintf("DataSize_%s_progressive_%s_%s.csv", file_kind, tag, start_tag))
  }
  split_t <- switch(split_method, random = "rnd", similarity = "simi", kennard_stone = "ks")
  start_t <- switch(
    initial_sampling,
    stratified_min = "strat",
    "rstart"
  )
  bal_t <- if (balance_method == "smote") "smote" else "classweight"
  kind_t <- switch(
    file_kind,
    metrics_raw = "raw",
    metrics_summary = "sum",
    best_models = "best",
    model_rankings = "rank",
    file_kind
  )
  sprintf("prog_%s_%s_%s_%s.csv", kind_t, split_t, bal_t, start_t)
}

win_long_path <- function(path) {
  if (.Platform$OS.type != "windows") {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  path <- normalizePath(path, winslash = "\\", mustWork = FALSE)
  if (grepl("^\\\\\\?\\\\", path)) {
    return(path)
  }
  if (nchar(path) < 248L) {
    return(path)
  }
  paste0("\\\\?\\", path)
}

fit_polr_weighted <- function(formula, data, weights) {
  if (length(weights) != nrow(data)) {
    stop("weights must have length nrow(data)", call. = FALSE)
  }
  formula_chr <- as.character(formula)[1]
  if (!nzchar(formula_chr)) return(NULL)
  f <- tryCatch(stats::as.formula(formula_chr), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  num.of.vars <- stringi::stri_count(formula_chr, fixed = "+")
  start <- c(rep(0, num.of.vars + 2), 1)
  success <- FALSE
  max_tries <- max(50L, num.of.vars + 25L)
  ntry <- 0L
  while (!success) {
    ntry <- ntry + 1L
    if (ntry > max_tries) break
    res <- tryCatch(
      list(
        ok = TRUE,
        model = MASS::polr(
          f,
          data = data,
          weights = weights,
          Hess = TRUE,
          start = start,
          control = list(maxit = 100)
        )
      ),
      error = function(e) list(ok = FALSE)
    )
    if (isTRUE(res$ok)) {
      success <- TRUE
      return(res$model)
    }
    start <- c(0, start)
  }
  NULL
}

make_class_weights_vec <- function(y) {
  y <- as.factor(y)
  tab <- table(y)
  k <- nlevels(y)
  n <- length(y)
  fac <- as.character(y)
  as.numeric(n / (k * as.numeric(tab[fac])))
}

apply_smote_train <- function(train_data, resp = "class") {
  if (!requireNamespace("recipes", quietly = TRUE) || !requireNamespace("themis", quietly = TRUE)) {
    stop("For balance_method='smote', install.packages(c('recipes', 'themis'))")
  }
  resp <- gsub("`", "", resp, fixed = TRUE)
  if (!resp %in% names(train_data)) stop("response column not found in train_data")
  pred_cols <- setdiff(names(train_data), resp)
  if (!length(pred_cols)) stop("SMOTE requires at least one predictor column besides the response")
  bq <- function(nm) paste0("`", gsub("`", "", nm, fixed = TRUE), "`")
  f <- stats::reformulate(
    termlabels = vapply(pred_cols, bq, character(1), USE.NAMES = FALSE),
    response = bq(resp)
  )
  rec0 <- recipes::recipe(f, data = train_data)
  rec1 <- themis::step_smote(rec0, recipes::all_outcomes(), over_ratio = 1)
  rec2 <- recipes::prep(rec1, training = train_data)
  recipes::bake(rec2, new_data = NULL)
}

prepare_modeling_train <- function(train_data, balance_method = BALANCE_METHODS) {
  balance_method <- match.arg(balance_method, BALANCE_METHODS)
  train_orig <- train_data
  if (balance_method == "smote") {
    train_for_fit <- tryCatch(
      apply_smote_train(train_orig),
      error = function(e) {
        warning("SMOTE failed; using unbalanced train. ", conditionMessage(e), call. = FALSE)
        train_orig
      }
    )
    train_for_fit$class <- factor(train_for_fit$class, levels = levels(train_orig$class))
    return(list(
      train_data = train_orig,
      train_for_fit = train_for_fit,
      case_weights = NULL,
      balance_method = balance_method
    ))
  }
  if (balance_method == "class_weight") {
    return(list(
      train_data = train_orig,
      train_for_fit = train_orig,
      case_weights = make_class_weights_vec(train_orig$class),
      balance_method = balance_method
    ))
  }
  list(
    train_data = train_orig,
    train_for_fit = train_orig,
    case_weights = NULL,
    balance_method = balance_method
  )
}

progressive_output_paths <- function(output_dir,
                                     split_method = SPLIT_METHODS,
                                     initial_sampling = INITIAL_SAMPLING_METHODS,
                                     balance_method = BALANCE_METHODS,
                                     inprogress = FALSE) {
  balance_method <- match.arg(balance_method, BALANCE_METHODS)
  split_method <- match.arg(split_method, SPLIT_METHODS)
  initial_sampling <- match.arg(initial_sampling, INITIAL_SAMPLING_METHODS)
  if (isTRUE(inprogress)) {
    bal_suffix <- switch(
      balance_method,
      none = "plain",
      smote = "smote",
      class_weight = "classweight"
    )
    return(list(
      raw = file.path(output_dir, sprintf("progressive_raw_%s_inprogress.csv", bal_suffix)),
      summary = file.path(output_dir, sprintf("progressive_summary_%s_inprogress.csv", bal_suffix)),
      best_models = file.path(
        output_dir, sprintf("progressive_best_models_%s_inprogress.csv", bal_suffix)
      ),
      model_rankings = file.path(
        output_dir, sprintf("progressive_rankings_%s_inprogress.csv", bal_suffix)
      )
    ))
  }
  list(
    raw = file.path(output_dir, progressive_csv_basename(
      "metrics_raw", split_method, initial_sampling, balance_method
    )),
    summary = file.path(output_dir, progressive_csv_basename(
      "metrics_summary", split_method, initial_sampling, balance_method
    )),
    best_models = file.path(output_dir, progressive_csv_basename(
      "best_models", split_method, initial_sampling, balance_method
    )),
    model_rankings = file.path(output_dir, progressive_csv_basename(
      "model_rankings", split_method, initial_sampling, balance_method
    ))
  )
}

finalize_one_data_size <- function(target_size,
                                   sk,
                                   top10_by_size,
                                   eval_by_size_rep,
                                   split_meta_by_size_rep,
                                   candidate_formulas_by_size,
                                   n_repeats) {
  top10_list <- top10_by_size[[sk]]
  if (is.null(top10_list) || length(top10_list) == 0) {
    return(NULL)
  }

  appearance_summary <- summarize_top10_appearances(top10_list)
  if (nrow(appearance_summary) == 0) {
    return(NULL)
  }

  best_formula <- normalize_formula_chr(appearance_summary$formula[1L])
  n_cand <- length(candidate_formulas_by_size[[sk]])
  if (is.null(n_cand) || is.na(n_cand)) n_cand <- 0L

  best_model <- data.frame(
    data_size = as.integer(target_size),
    best_formula = best_formula,
    top10_count = appearance_summary$top10_count[1L],
    n_candidates = as.integer(n_cand),
    stringsAsFactors = FALSE
  )

  appearance_summary$data_size <- as.integer(target_size)
  model_rankings <- appearance_summary[, c("data_size", "formula", "top10_count"), drop = FALSE]

  rep_evals <- eval_by_size_rep[[sk]]
  rep_meta <- split_meta_by_size_rep[[sk]]
  metrics_raw_rows <- vector("list", n_repeats)

  for (rep_i in seq_len(n_repeats)) {
    eval_df <- rep_evals[[rep_i]]
    meta <- rep_meta[[rep_i]]
    n_tr <- if (is.null(meta)) NA_integer_ else meta$n_train
    n_te <- if (is.null(meta)) NA_integer_ else meta$n_test

    if (is.null(eval_df) || nrow(eval_df) == 0) {
      metrics_raw_rows[[rep_i]] <- data.frame(
        data_size = as.integer(target_size),
        repeat_id = rep_i,
        best_formula = best_formula,
        n_train = n_tr,
        n_test = n_te,
        train_accuracy = NA_real_,
        test_accuracy = NA_real_,
        external_accuracy = NA_real_,
        smallest_group_recall = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    hit <- eval_df[normalize_formula_chr(eval_df$formula) == best_formula, , drop = FALSE]
    metrics_raw_rows[[rep_i]] <- data.frame(
      data_size = as.integer(target_size),
      repeat_id = rep_i,
      best_formula = best_formula,
      n_train = n_tr,
      n_test = n_te,
      train_accuracy = if (nrow(hit) == 0) NA_real_ else hit$train_accuracy[1],
      test_accuracy = if (nrow(hit) == 0) NA_real_ else hit$test_accuracy[1],
      external_accuracy = if (nrow(hit) == 0) NA_real_ else hit$external_accuracy[1],
      smallest_group_recall = if (nrow(hit) == 0 || !"smallest_group_recall" %in% names(hit)) {
        NA_real_
      } else {
        hit$smallest_group_recall[1]
      },
      stringsAsFactors = FALSE
    )
  }

  metrics_raw <- do.call(rbind, metrics_raw_rows)
  rownames(metrics_raw) <- NULL

  metrics_summary <- data.frame(
    data_size = as.integer(target_size),
    best_formula = best_formula,
    train_accuracy = safe_mean(metrics_raw$train_accuracy),
    train_accuracy_sd = safe_sd(metrics_raw$train_accuracy),
    test_accuracy = safe_mean(metrics_raw$test_accuracy),
    test_accuracy_sd = safe_sd(metrics_raw$test_accuracy),
    external_accuracy = safe_mean(metrics_raw$external_accuracy),
    external_accuracy_sd = safe_sd(metrics_raw$external_accuracy),
    smallest_group_recall = safe_mean(metrics_raw$smallest_group_recall),
    smallest_group_recall_sd = safe_sd(metrics_raw$smallest_group_recall),
    stringsAsFactors = FALSE
  )

  list(
    model_rankings = model_rankings,
    best_model = best_model,
    metrics_raw = metrics_raw,
    metrics_summary = metrics_summary
  )
}

summarize_metrics_from_raw <- function(metrics_raw, target_sizes, best_models = NULL) {
  sizes <- as.integer(target_sizes)
  has_recall <- "smallest_group_recall" %in% names(metrics_raw)
  out <- data.frame(
    data_size = sizes,
    train_accuracy = vapply(sizes, function(s) {
      safe_mean(metrics_raw$train_accuracy[metrics_raw$data_size == s])
    }, numeric(1)),
    train_accuracy_sd = vapply(sizes, function(s) {
      safe_sd(metrics_raw$train_accuracy[metrics_raw$data_size == s])
    }, numeric(1)),
    test_accuracy = vapply(sizes, function(s) {
      safe_mean(metrics_raw$test_accuracy[metrics_raw$data_size == s])
    }, numeric(1)),
    test_accuracy_sd = vapply(sizes, function(s) {
      safe_sd(metrics_raw$test_accuracy[metrics_raw$data_size == s])
    }, numeric(1)),
    external_accuracy = vapply(sizes, function(s) {
      safe_mean(metrics_raw$external_accuracy[metrics_raw$data_size == s])
    }, numeric(1)),
    external_accuracy_sd = vapply(sizes, function(s) {
      safe_sd(metrics_raw$external_accuracy[metrics_raw$data_size == s])
    }, numeric(1)),
    smallest_group_recall = if (has_recall) {
      vapply(sizes, function(s) {
        safe_mean(metrics_raw$smallest_group_recall[metrics_raw$data_size == s])
      }, numeric(1))
    } else {
      rep(NA_real_, length(sizes))
    },
    smallest_group_recall_sd = if (has_recall) {
      vapply(sizes, function(s) {
        safe_sd(metrics_raw$smallest_group_recall[metrics_raw$data_size == s])
      }, numeric(1))
    } else {
      rep(NA_real_, length(sizes))
    },
    stringsAsFactors = FALSE
  )
  if (!is.null(best_models) && nrow(best_models) > 0) {
    best_by_size <- stats::setNames(best_models$best_formula, best_models$data_size)
    out$best_formula <- vapply(out$data_size, function(s) {
      bf <- best_by_size[[as.character(s)]]
      if (is.null(bf)) NA_character_ else bf
    }, character(1))
    out <- out[, c(
      "data_size", "best_formula",
      "train_accuracy", "train_accuracy_sd",
      "test_accuracy", "test_accuracy_sd",
      "external_accuracy", "external_accuracy_sd",
      "smallest_group_recall", "smallest_group_recall_sd"
    )]
  }
  out
}

ensure_dir_exists <- function(path, what = "output directory") {
  path <- trimws(as.character(path)[1L])
  if (!nzchar(path) || is.na(path)) {
    stop("Invalid ", what, ": empty or NA path.", call. = FALSE)
  }
  if (!dir.exists(path)) {
    ok <- suppressWarnings(dir.create(path, recursive = TRUE))
    if (!ok && !dir.exists(path)) {
      stop(
        "Could not create ", what, ":\n  ", path,
        "\n(Check that parent folders exist and you have write permission.)",
        call. = FALSE
      )
    }
  }
  if (!dir.exists(path)) {
    stop(what, " does not exist:\n  ", path, call. = FALSE)
  }
  tryCatch(
    normalizePath(path, winslash = "/", mustWork = TRUE),
    error = function(e) normalizePath(path, winslash = "/", mustWork = FALSE)
  )
}

write_df_csv <- function(df, path) {
  if (is.null(df)) {
    stop("Cannot write NULL data frame to: ", path, call. = FALSE)
  }
  path <- as.character(path)[1L]
  out_dir <- dirname(path)
  if (!nzchar(out_dir) || identical(out_dir, ".") || is.na(out_dir)) {
    stop("Invalid CSV path (missing parent directory): ", path, call. = FALSE)
  }
  out_dir <- ensure_dir_exists(out_dir, what = "directory for progressive CSV exports")
  path <- file.path(out_dir, basename(path))
  write_path <- win_long_path(path)

  write_once <- function(target) {
    utils::write.csv(df, file = target, row.names = FALSE)
  }

  err <- NULL
  ok <- tryCatch({
    write_once(write_path)
    TRUE
  }, error = function(e) {
    err <<- conditionMessage(e)
    FALSE
  })
  exists_at <- function(p) file.exists(p) || file.exists(path)
  if (!ok || !exists_at(write_path)) {
    tmp <- file.path(out_dir, paste0("_tmp_", basename(path)))
    tmp_long <- win_long_path(tmp)
    ok2 <- tryCatch({
      write_once(tmp_long)
      if (file.exists(path)) unlink(path)
      if (file.exists(write_path) && write_path != path) {
        try(unlink(write_path), silent = TRUE)
      }
      renamed <- file.rename(tmp, path)
      if (!renamed) {
        file.copy(tmp, path, overwrite = TRUE)
        unlink(tmp)
      }
      TRUE
    }, error = function(e2) {
      err <<- paste(err, "| fallback:", conditionMessage(e2))
      FALSE
    })
    if (!ok2 || !file.exists(path)) {
      stop(
        "Failed to write CSV.\n  Target: ", path,
        "\n  Path length: ", nchar(path),
        "\n  Directory exists: ", dir.exists(out_dir),
        "\n  Error: ", err,
        call. = FALSE
      )
    }
  }
  invisible(path)
}

persist_progressive_metrics <- function(paths,
                                        metrics_raw,
                                        metrics_summary,
                                        best_models,
                                        model_rankings) {
  objs <- list(
    raw = metrics_raw,
    summary = metrics_summary,
    best_models = best_models,
    model_rankings = model_rankings
  )
  for (label in names(objs)) {
    tryCatch(
      write_df_csv(objs[[label]], paths[[label]]),
      error = function(e) {
        stop("Failed writing ", label, " to ", paths[[label]], ": ", conditionMessage(e), call. = FALSE)
      }
    )
    if (!file.exists(paths[[label]])) {
      stop("CSV not found after write: ", paths[[label]], call. = FALSE)
    }
  }
  invisible(paths)
}

write_progressive_inprogress_snapshot <- function(inprogress_paths,
                                                  metrics_raw,
                                                  metrics_summary,
                                                  best_models,
                                                  model_rankings,
                                                  through_label = "run") {
  persist_progressive_metrics(
    paths = inprogress_paths,
    metrics_raw = metrics_raw,
    metrics_summary = metrics_summary,
    best_models = best_models,
    model_rankings = model_rankings
  )
  message(sprintf(
    paste0(
      "In-progress CSVs updated (%s):\n",
      "  summary: %s\n  raw: %s\n  best_models: %s\n  rankings: %s"
    ),
    through_label,
    inprogress_paths$summary,
    inprogress_paths$raw,
    inprogress_paths$best_models,
    inprogress_paths$model_rankings
  ))
  flush.console()
  invisible(inprogress_paths)
}

print_size_step_results <- function(target_size, step_res) {
  message(sprintf("\n========== Completed data_size=%d (%d repeats) ==========", target_size, nrow(step_res$metrics_raw)))
  message("--- Model rankings (top-10 appearances) ---")
  print(step_res$model_rankings)
  message("--- Best model (highest top-10 count) ---")
  print(step_res$best_model)
  message("--- Mean metrics for best model across repeats ---")
  print(step_res$metrics_summary)
  message("====================================================\n")
}

submodel_bounds_for_size <- function(data_size) {
  sz <- as.integer(data_size)[1L]
  if (is.na(sz)) sz <- 25L
  if (sz <= 25L) {
    list(min = 3L, max = 3L)
  } else {
    list(min = 4L, max = 4L)
  }
}

run_submodel_search <- function(train_data,
                                data_size,
                                ordinal = TRUE,
                                submodel_min = NULL,
                                submodel_max = NULL,
                                n_candidates = 200L,
                                case_weights = NULL) {
  if (is.null(submodel_min) || is.null(submodel_max)) {
    bounds <- submodel_bounds_for_size(data_size)
    submodel_min <- bounds$min
    submodel_max <- bounds$max
  }
  if (nrow(train_data) < 3 || length(unique(as.character(train_data$class))) < 2) {
    return(NULL)
  }
  out <- tryCatch(
    sub_model_log(
      data = train_data,
      min = submodel_min,
      max = submodel_max,
      ordinal = ordinal,
      weights = case_weights,
      top_n = n_candidates
    ),
    error = function(e) NULL
  )
  if (!is.null(out) && nrow(out) < n_candidates) {
    message(sprintf(
      "sub_model_log returned %d models (requested %d). Reload local rxn.cond.class (top_n support) or lower n_candidates.",
      nrow(out), n_candidates
    ))
  }
  out
}

build_external_eval_data <- function(all_data, selected_idx, external_data = NULL) {
  align_to_template <- function(df, template_names, class_levels = NULL) {
    if (is.null(df)) return(NULL)
    if (nrow(df) == 0) {
      out <- as.data.frame(matrix(nrow = 0, ncol = length(template_names)))
      colnames(out) <- template_names
      if ("class" %in% template_names && !is.null(class_levels)) {
        out$class <- factor(character(0), levels = class_levels)
      }
      return(out)
    }

    missing_cols <- setdiff(template_names, colnames(df))
    if (length(missing_cols) > 0) {
      for (nm in missing_cols) {
        df[[nm]] <- NA
      }
    }
    df <- df[, template_names, drop = FALSE]
    if ("class" %in% template_names && !is.null(class_levels)) {
      df$class <- factor(as.character(df$class), levels = class_levels)
    }
    df
  }

  unsampled_data <- all_data[setdiff(seq_len(nrow(all_data)), selected_idx), , drop = FALSE]
  template_names <- colnames(all_data)
  class_levels <- if ("class" %in% template_names) levels(all_data$class) else NULL
  unsampled_data <- align_to_template(unsampled_data, template_names, class_levels)

  if (is.null(external_data) || nrow(external_data) == 0) {
    return(unsampled_data)
  }
  external_data <- align_to_template(external_data, template_names, class_levels)
  if (nrow(unsampled_data) == 0) {
    return(external_data)
  }

  eval_data <- rbind(external_data, unsampled_data)
  row.names(eval_data) <- make.unique(as.character(row.names(eval_data)))
  eval_data
}

sample_subset_indices <- function(data,
                                  subset_size,
                                  seed = NULL,
                                  method = c("random", "stratified_min"),
                                  min_per_class = 5) {
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)
  n_total <- nrow(data)
  if (n_total == 0) return(integer(0))

  subset_size <- as.integer(subset_size)[1]
  if (is.na(subset_size) || subset_size <= 0) return(integer(0))
  subset_size <- min(subset_size, n_total)

  if (method == "random") {
    return(sort(sample(seq_len(n_total), size = subset_size, replace = FALSE)))
  }

  selected <- integer(0)
  cls <- levels(data$class)
  for (cl in cls) {
    idx_cl <- which(as.character(data$class) == cl)
    if (length(idx_cl) == 0) next
    take_cl <- min(as.integer(min_per_class), length(idx_cl), subset_size - length(selected))
    if (take_cl > 0L) {
      selected <- c(selected, sample(idx_cl, size = take_cl, replace = FALSE))
    }
    if (length(selected) >= subset_size) break
  }

  selected <- unique(selected)
  if (length(selected) < subset_size) {
    remaining <- setdiff(seq_len(n_total), selected)
    extra_n <- min(subset_size - length(selected), length(remaining))
    if (extra_n > 0L) {
      selected <- c(selected, sample(remaining, size = extra_n, replace = FALSE))
    }
  }
  sort(unique(selected))
}

sample_indices_for_target_size <- function(data,
                                           target_size,
                                           seed = NULL,
                                           subset_sampling = c("random", "stratified_min"),
                                           min_per_class = 5) {
  subset_sampling <- match.arg(subset_sampling)
  sample_subset_indices(
    data = data,
    subset_size = as.integer(target_size)[1L],
    seed = seed,
    method = subset_sampling,
    min_per_class = min_per_class
  )
}

# Progressive data-size study: sample subsets at each size, search formulas (repeat 1),
# re-evaluate locked candidates on new samples (repeats 2..n), pick stable best model per size.
#
# Arguments:
#   data              Training pool (class + predictors)
#   external_data     Optional held-out set for external_accuracy
#   split_method      Train/holdout split: random, similarity, or kennard_stone
#   initial_sampling  How each subset is drawn: random or stratified_min
#   min_per_class     Minimum rows per class in stratified_min sampling
#   start_size        First target subset size
#   max_size          Last target subset size
#   step_add          Subset size increment between steps
#   train_p           Holdout fraction for random split (when split_method=random)
#   random_train_p    Holdout fraction override for similarity/KS splits (default: train_p)
#   target_train_n    Fixed train size cap for similarity/KS splits (default: full sampled set)
#   n_repeats         Independent repeats per target size
#   ordinal           TRUE = polr, FALSE = multinom
#   n_candidates      Formulas from sub_model_log on repeat 1
#   top_k_per_repeat  Top formulas kept per repeat for stability ranking
#   balance_method    none, smote, or class_weight (training split only)
#   output_dir        Folder for in-progress CSV snapshots
#   write_inprogress_csv  Save partial results after each size
#   seed              Base random seed
evaluate_progressive_data_size <- function(data,
                                           external_data = NULL,
                                           split_method = c("random", "similarity", "kennard_stone"),
                                           initial_sampling = c("random", "stratified_min"),
                                           min_per_class = 5,
                                           start_size = 20,
                                           max_size = 55,
                                           step_add = 5,
                                           train_p = 0.75,
                                           random_train_p = NULL,
                                           target_train_n = NULL,
                                           n_repeats = 100,
                                           ordinal = TRUE,
                                           n_candidates = 200L,
                                           top_k_per_repeat = 10L,
                                           balance_method = c("none", "smote", "class_weight"),
                                           output_dir = NULL,
                                           write_inprogress_csv = TRUE,
                                           seed = 1) {
  split_method <- match.arg(split_method)
  initial_sampling <- match.arg(initial_sampling)
  balance_method <- match.arg(balance_method)
  if (is.null(random_train_p)) random_train_p <- train_p
  if (split_method == "kennard_stone") {
    ensure_kennard_stone_helpers(re_source = TRUE)
  }
  if (!is.factor(data$class)) data$class <- factor(data$class)
  n_total <- nrow(data)
  class_counts <- table(data$class)
  smallest_class <- names(class_counts)[which.min(as.integer(class_counts))][1]
  target_sizes <- seq(start_size, max_size, by = step_add)

  size_key <- function(sz) as.character(as.integer(sz))
  candidate_formulas_by_size <- list()
  eval_by_size_rep <- list()
  split_meta_by_size_rep <- list()
  top10_by_size <- setNames(
    lapply(target_sizes, function(sz) vector("list", n_repeats)),
    vapply(target_sizes, size_key, character(1))
  )

  completed_sizes <- integer(0)
  best_models_accum <- list()
  model_rankings_accum <- list()
  metrics_raw_accum <- list()

  inprogress_paths <- NULL
  if (write_inprogress_csv && !is.null(output_dir) && nzchar(output_dir)) {
    output_dir <- ensure_dir_exists(output_dir, what = "output_dir (in-progress CSVs)")
    inprogress_paths <- progressive_output_paths(
      output_dir = output_dir,
      split_method = split_method,
      initial_sampling = initial_sampling,
      balance_method = balance_method,
      inprogress = TRUE
    )
    bal_lab <- switch(balance_method, none = "plain", smote = "smote", class_weight = "classweight")
    message(sprintf(
      "In-progress saves enabled (balance=%s; files end with _%s_inprogress.csv):",
      balance_method, bal_lab
    ))
    message(sprintf("  %s", inprogress_paths$summary))
    message(sprintf("  %s", inprogress_paths$raw))
    message(sprintf("  %s", inprogress_paths$best_models))
    message(sprintf("  %s", inprogress_paths$model_rankings))
  }

  message(sprintf(
    "Balance method: %s (SMOTE/class weights apply to train split only; metrics on real holdout/external/full data)",
    balance_method
  ))
  message(sprintf(
    "Split: %s | Initial subset sampling: %s (min_per_class=%d)",
    split_method, initial_sampling, as.integer(min_per_class)
  ))
  if (split_method == "kennard_stone") {
    message("Kennard–Stone split: features scaled globally; 'flag' excluded from distance.")
  }

  total_expected_steps <- as.integer(length(target_sizes) * n_repeats)
  global_step <- 0L

  for (step_in_sizes in seq_along(target_sizes)) {
    target_size <- as.integer(target_sizes[step_in_sizes])
    sk <- size_key(target_size)
    candidate_formulas_by_size[[sk]] <- NULL
    eval_by_size_rep[[sk]] <- vector("list", n_repeats)
    split_meta_by_size_rep[[sk]] <- vector("list", n_repeats)

    message(sprintf(
      "[Data size %d/%d] target_size=%d — independent sampling across %d repeats",
      step_in_sizes, length(target_sizes), target_size, n_repeats
    ))

    for (rep_i in seq_len(n_repeats)) {
      global_step <- global_step + 1L
      sample_seed <- seed + target_size * 10000L + rep_i
      selected_idx <- sample_indices_for_target_size(
        data = data,
        target_size = target_size,
        seed = sample_seed,
        subset_sampling = initial_sampling,
        min_per_class = min_per_class
      )

      if (length(selected_idx) < target_size) {
        warning(sprintf(
          "data_size=%d repeat=%d skipped: sampled %d rows (need %d).",
          target_size, rep_i, length(selected_idx), target_size
        ))
        eval_by_size_rep[[sk]][[rep_i]] <- NULL
        top10_by_size[[sk]][[rep_i]] <- character(0)
        next
      }

      sampled_data <- data[selected_idx, , drop = FALSE]
      sampled_data$class <- factor(sampled_data$class, levels = levels(data$class))
      n_sampled <- nrow(sampled_data)
      left_out_idx <- setdiff(seq_len(n_total), selected_idx)

      message(sprintf(
        "[Step %d/%d] data_size=%d repeat=%d/%d n_sampled=%d (unsampled corpus rows=%d)",
        global_step, total_expected_steps, target_size, rep_i, n_repeats, n_sampled, length(left_out_idx)
      ))

      split_seed <- seed + target_size * 1000L + rep_i
      split_obj <- if (split_method == "similarity") {
        base_n <- if (is.null(target_train_n)) nrow(sampled_data) else as.integer(target_train_n)[1]
        if (is.na(base_n)) base_n <- nrow(sampled_data)
        target_n <- as.integer(round(base_n * random_train_p))
        target_n <- max(1L, min(target_n, nrow(sampled_data) - 1L))
        split_with_similarity_sampler(sampled_data, seed = split_seed, target_train_n = target_n)
      } else if (split_method == "kennard_stone") {
        base_n <- if (is.null(target_train_n)) nrow(sampled_data) else as.integer(target_train_n)[1]
        if (is.na(base_n)) base_n <- nrow(sampled_data)
        target_n <- as.integer(round(base_n * random_train_p))
        target_n <- max(1L, min(target_n, nrow(sampled_data) - 1L))
        split_with_kennard_stone_sampler(
          sampled_data,
          seed = split_seed,
          target_train_n = target_n,
          min_per_class = min_per_class
        )
      } else {
        split_with_random_partition(sampled_data, seed = split_seed, train_p = random_train_p)
      }

      class_levels <- levels(data$class)
      train_data <- align_class_levels(split_obj$train_data, class_levels)
      holdout_test_data <- align_class_levels(split_obj$test_data, class_levels)
      holdout_test_data$class <- factor(
        as.character(holdout_test_data$class),
        levels = class_levels
      )

      split_meta_by_size_rep[[sk]][[rep_i]] <- list(
        n_train = nrow(train_data),
        n_test = nrow(holdout_test_data)
      )

      message(sprintf(
        "  Split sizes: train_n=%d, holdout_test_n=%d (test accuracy uses holdout only; external uses external+unsampled)",
        nrow(train_data), nrow(holdout_test_data)
      ))

      if (nrow(train_data) < 3 || length(unique(as.character(train_data$class))) < 2) {
        eval_by_size_rep[[sk]][[rep_i]] <- NULL
        top10_by_size[[sk]][[rep_i]] <- character(0)
        message(sprintf(
          "  Insufficient train diversity (n_train=%d, classes=%d). Skipped.",
          nrow(train_data), length(unique(as.character(train_data$class)))
        ))
        next
      }

      train_prep <- prepare_modeling_train(train_data, balance_method = balance_method)
      if (balance_method == "smote" && nrow(train_prep$train_for_fit) != nrow(train_prep$train_data)) {
        message(sprintf(
          "  SMOTE train rows: %d (original train n=%d)",
          nrow(train_prep$train_for_fit), nrow(train_prep$train_data)
        ))
      }

      mcfadden_lookup <- NULL
      if (rep_i == 1L) {
        sm_bounds <- submodel_bounds_for_size(target_size)
        message(sprintf(
          "  Repeat 1: sub_model_log (top %d, %d predictors) on sampled n=%d, train n=%d (%s)...",
          n_candidates, sm_bounds$min, target_size, nrow(train_prep$train_for_fit), balance_method
        ))
        search_models <- run_submodel_search(
          train_data = train_prep$train_for_fit,
          data_size = target_size,
          ordinal = ordinal,
          n_candidates = n_candidates,
          case_weights = train_prep$case_weights
        )
        if (is.null(search_models) || nrow(search_models) == 0) {
          warning(sprintf("sub_model_log returned no models at data_size=%d.", target_size))
          eval_by_size_rep[[sk]][[rep_i]] <- NULL
          top10_by_size[[sk]][[rep_i]] <- character(0)
          next
        }
        mcf_name <- mcfadden_col_name(search_models)
        candidate_formulas_by_size[[sk]] <- normalize_formula_chr(search_models$formula)
        if (!is.na(mcf_name)) {
          mcfadden_lookup <- stats::setNames(
            as.numeric(search_models[[mcf_name]]),
            candidate_formulas_by_size[[sk]]
          )
        }
        message(sprintf(
          "  Locked %d candidate formulas for data_size=%d (used in repeats 2-%d).",
          length(candidate_formulas_by_size[[sk]]), target_size, n_repeats
        ))
      } else {
        if (is.null(candidate_formulas_by_size[[sk]]) || length(candidate_formulas_by_size[[sk]]) == 0) {
          warning(sprintf(
            "No candidates from repeat 1 at data_size=%d; repeat %d skipped.",
            target_size, rep_i
          ))
          next
        }
        message(sprintf(
          "  Repeat %d: evaluating %d locked candidates...",
          rep_i, length(candidate_formulas_by_size[[sk]])
        ))
      }

      eval_df <- evaluate_formula_candidates(
        formulas = candidate_formulas_by_size[[sk]],
        train_data = train_prep$train_data,
        holdout_test_data = holdout_test_data,
        all_data = data,
        sampled_idx = selected_idx,
        external_data = external_data,
        smallest_class = smallest_class,
        train_for_fit = train_prep$train_for_fit,
        case_weights = train_prep$case_weights,
        mcfadden_train = mcfadden_lookup,
        ordinal = ordinal
      )
      eval_df$formula <- normalize_formula_chr(eval_df$formula)
      eval_by_size_rep[[sk]][[rep_i]] <- eval_df
      top10_formulas <- rank_top_formulas(eval_df, top_k = top_k_per_repeat)
      top10_by_size[[sk]][[rep_i]] <- top10_formulas

      n_ok <- sum(!is.na(eval_df$test_accuracy))
      message(sprintf(
        "  Repeat %d done: %d/%d candidates fitted; top-%d selected for stability count.",
        rep_i, n_ok, nrow(eval_df), top_k_per_repeat
      ))
      if (balance_method == "class_weight" && n_ok == 0L && nrow(eval_df) > 0L) {
        warning(
          sprintf(
            paste0(
              "data_size=%d repeat=%d: all %d candidate fits failed (class_weight). ",
              "Re-source this script after devtools::load_all(rxn.cond.class); ",
              "compare with Model Construction_SMOTE_Class_Weighting.R on the same split."
            ),
            target_size, rep_i, nrow(eval_df)
          ),
          call. = FALSE
        )
      }

      if (rep_i > 1L && nrow(eval_df) > 0L) {
        ex <- eval_df[1L, , drop = FALSE]
        mcf_lab <- if (is.na(ex$mcfadden_train[1])) "NA" else sprintf("%.3f", ex$mcfadden_train[1])
        message(sprintf(
          paste0(
            "  Repeat %d — locked candidate #1: train=%.3f, test=%.3f, external+unsampled=%.3f, recall_smallest=%.3f, McFadden=%s",
            " | formula: %s"
          ),
          rep_i,
          ex$train_accuracy[1], ex$test_accuracy[1], ex$external_accuracy[1],
          ex$smallest_group_recall[1],
          mcf_lab,
          ex$formula[1]
        ))
      }
    }

    step_res <- finalize_one_data_size(
      target_size = target_size,
      sk = sk,
      top10_by_size = top10_by_size,
      eval_by_size_rep = eval_by_size_rep,
      split_meta_by_size_rep = split_meta_by_size_rep,
      candidate_formulas_by_size = candidate_formulas_by_size,
      n_repeats = n_repeats
    )
    if (is.null(step_res)) {
      warning(sprintf("No stability results to aggregate for data_size=%d.", target_size))
      next
    }

    print_size_step_results(target_size, step_res)
    completed_sizes <- c(completed_sizes, target_size)
    best_models_accum[[length(best_models_accum) + 1L]] <- step_res$best_model
    model_rankings_accum[[length(model_rankings_accum) + 1L]] <- step_res$model_rankings
    metrics_raw_accum[[length(metrics_raw_accum) + 1L]] <- step_res$metrics_raw

    if (!is.null(inprogress_paths)) {
      best_models <- do.call(rbind, best_models_accum)
      model_rankings <- do.call(rbind, model_rankings_accum)
      metrics_raw <- do.call(rbind, metrics_raw_accum)
      metrics_summary <- summarize_metrics_from_raw(metrics_raw, completed_sizes, best_models)
      write_progressive_inprogress_snapshot(
        inprogress_paths = inprogress_paths,
        metrics_raw = metrics_raw,
        metrics_summary = metrics_summary,
        best_models = best_models,
        model_rankings = model_rankings,
        through_label = sprintf("through data_size=%d", target_size)
      )
    }
  }

  if (length(best_models_accum) > 0) {
    best_models <- do.call(rbind, best_models_accum)
    rownames(best_models) <- NULL
  } else {
    best_models <- data.frame(
      data_size = integer(0),
      best_formula = character(0),
      top10_count = integer(0),
      n_candidates = integer(0),
      stringsAsFactors = FALSE
    )
  }

  if (length(model_rankings_accum) > 0) {
    model_rankings <- do.call(rbind, model_rankings_accum)
    rownames(model_rankings) <- NULL
  } else {
    model_rankings <- data.frame(
      data_size = integer(0),
      formula = character(0),
      top10_count = integer(0),
      stringsAsFactors = FALSE
    )
  }

  if (length(metrics_raw_accum) > 0) {
    metrics_raw <- do.call(rbind, metrics_raw_accum)
    rownames(metrics_raw) <- NULL
  } else {
    metrics_raw <- data.frame(
      data_size = integer(0),
      repeat_id = integer(0),
      best_formula = character(0),
      n_train = integer(0),
      n_test = integer(0),
      train_accuracy = numeric(0),
      test_accuracy = numeric(0),
      external_accuracy = numeric(0),
      smallest_group_recall = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  sizes_for_summary <- if (length(completed_sizes) > 0) completed_sizes else target_sizes
  metrics_summary <- summarize_metrics_from_raw(metrics_raw, sizes_for_summary, best_models)

  if (length(completed_sizes) == 0L) {
    warning(
      "No data sizes produced stability results (top-10 empty every time). ",
      "Check sub_model_log / class_weight fits. In-progress CSVs were not written.",
      call. = FALSE
    )
  } else if (!is.null(inprogress_paths)) {
    write_progressive_inprogress_snapshot(
      inprogress_paths = inprogress_paths,
      metrics_raw = metrics_raw,
      metrics_summary = metrics_summary,
      best_models = best_models,
      model_rankings = model_rankings,
      through_label = "end of evaluate_progressive_data_size"
    )
  }

  plot_df <- rbind(
    data.frame(data_size = metrics_summary$data_size, metric = "train_accuracy", value = metrics_summary$train_accuracy),
    data.frame(data_size = metrics_summary$data_size, metric = "test_accuracy", value = metrics_summary$test_accuracy),
    data.frame(data_size = metrics_summary$data_size, metric = "external_accuracy", value = metrics_summary$external_accuracy)
  )
  if ("smallest_group_recall" %in% names(metrics_summary)) {
    plot_df <- rbind(
      plot_df,
      data.frame(
        data_size = metrics_summary$data_size,
        metric = "smallest_group_recall",
        value = metrics_summary$smallest_group_recall
      )
    )
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes_string(x = "data_size", y = "value", color = "metric")) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, by = 0.1)) +
    ggplot2::labs(
      title = sprintf("Data-size evaluation — best model by stability (%s split)", split_method),
      x = "Data size",
      y = "Accuracy",
      color = "Metric"
    ) +
    ggplot2::theme_minimal()

  list(
    metrics_raw = metrics_raw,
    metrics_summary = metrics_summary,
    best_models = best_models,
    model_rankings = model_rankings,
    plot = p,
    balance_method = balance_method
  )
}

finalize_progressive_inprogress_to_final <- function(output_dir,
                                                     split_method = SPLIT_METHODS,
                                                     initial_sampling = INITIAL_SAMPLING_METHODS,
                                                     balance_method = BALANCE_METHODS) {
  output_dir <- ensure_dir_exists(output_dir, what = "output_dir (finalize in-progress)")
  paths_ip <- progressive_output_paths(
    output_dir, split_method, initial_sampling, balance_method,
    inprogress = TRUE
  )
  paths_fin <- progressive_output_paths(
    output_dir, split_method, initial_sampling, balance_method,
    inprogress = FALSE
  )
  for (label in names(paths_ip)) {
    ip <- paths_ip[[label]]
    fin <- paths_fin[[label]]
    if (!file.exists(ip)) {
      warning("Missing in-progress file, skipped: ", ip, call. = FALSE)
      next
    }
    df <- utils::read.csv(ip, stringsAsFactors = FALSE, check.names = FALSE)
    write_df_csv(df, fin)
    message(sprintf("Finalized %s -> %s", basename(ip), basename(fin)))
  }
  invisible(paths_fin)
}

save_progressive_outputs <- function(res,
                                     output_dir,
                                     split_method = SPLIT_METHODS,
                                     initial_sampling = INITIAL_SAMPLING_METHODS,
                                     balance_method = BALANCE_METHODS) {
  split_method <- match.arg(split_method, SPLIT_METHODS)
  initial_sampling <- match.arg(initial_sampling, INITIAL_SAMPLING_METHODS)
  balance_method <- match.arg(balance_method, BALANCE_METHODS)
  tag <- switch(split_method, random = "random", similarity = "simi", kennard_stone = "ks")
  start_tag <- switch(
    initial_sampling,
    stratified_min = "stratstart",
    "randomstart"
  )
  output_dir <- ensure_dir_exists(output_dir, what = "output_dir (final CSVs)")
  paths <- progressive_output_paths(
    output_dir = output_dir,
    split_method = split_method,
    initial_sampling = initial_sampling,
    balance_method = balance_method,
    inprogress = FALSE
  )
  raw_csv <- paths$raw
  summary_csv <- paths$summary
  best_models_csv <- paths$best_models
  model_rankings_csv <- paths$model_rankings
  plot_bal <- if (balance_method == "none") {
    ""
  } else if (balance_method == "smote") {
    "_smote"
  } else {
    "_classweight"
  }
  plot_png <- file.path(
    output_dir,
    sprintf("Progressive_accuracy_plot_%s%s_%s.png", tag, plot_bal, start_tag)
  )

  if (nrow(res$metrics_raw) == 0L) {
    warning(
      "metrics_raw is empty — the progressive run did not aggregate any size. ",
      "Final CSVs will contain headers/NA only; fix modeling before trusting results.",
      call. = FALSE
    )
  }

  write_df_csv(res$metrics_raw, raw_csv)
  write_df_csv(res$metrics_summary, summary_csv)
  write_df_csv(res$best_models, best_models_csv)
  write_df_csv(res$model_rankings, model_rankings_csv)
  ggplot2::ggsave(filename = plot_png, plot = res$plot, width = 8, height = 5, dpi = 300)

  inprogress_paths <- progressive_output_paths(
    output_dir = output_dir,
    split_method = split_method,
    initial_sampling = initial_sampling,
    balance_method = balance_method,
    inprogress = TRUE
  )
  if (nrow(res$metrics_raw) > 0L) {
    write_progressive_inprogress_snapshot(
      inprogress_paths = inprogress_paths,
      metrics_raw = res$metrics_raw,
      metrics_summary = res$metrics_summary,
      best_models = res$best_models,
      model_rankings = res$model_rankings,
      through_label = "mirrored from final save_progressive_outputs"
    )
  }

  message(sprintf("Saved raw metrics: %s", raw_csv))
  message(sprintf("Saved summary metrics: %s", summary_csv))
  message(sprintf("Saved best models: %s", best_models_csv))
  message(sprintf("Saved model rankings: %s", model_rankings_csv))
  message(sprintf("Saved plot: %s", plot_png))
  if (balance_method != "none") {
    message("Note: balanced runs use shorter output names.")
  }
}

env_csv_override <- function(var_names) {
  for (nm in var_names) {
    val <- Sys.getenv(nm, unset = "")
    if (nzchar(val) && file.exists(val)) {
      return(normalizePath(val, winslash = "/", mustWork = TRUE))
    }
  }
  NULL
}

load_training_data <- function(path) {
  for (fn in c("prepare_training_data", "prepare_deuteration_training_data")) {
    if (exists(fn, mode = "function")) {
      return(get(fn)(path))
    }
  }
  data <- data.frame(data.table::fread(path, check.names = FALSE))
  if (ncol(data) < 3L) {
    stop("expected tag/class columns in training data", call. = FALSE)
  }
  row.names(data) <- data[, 2L]
  data$class <- as.factor(data$class)
  data <- data[, -c(1L, 2L), drop = FALSE]
  pred_idx <- setdiff(seq_len(ncol(data)), which(names(data) == "class"))
  data[pred_idx] <- lapply(data[pred_idx], function(x) suppressWarnings(as.numeric(x)))
  data
}

find_training_csv <- function(start_dirs, max_up = 8L) {
  hit <- env_csv_override(c("TRAINING_CSV", "DEUTERATION_TRAINING_CSV"))
  if (!is.null(hit)) return(hit)

  search_from <- unique(normalizePath(start_dirs, winslash = "/", mustWork = FALSE))
  for (root in search_from) {
    cur <- root
    for (i in 0:max_up) {
      candidate <- file.path(cur, "Data", "Training_Data.csv")
      if (file.exists(candidate)) {
        return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  NULL
}

load_external_validation_aligned <- function(external_csv, ref_data) {
  if (is.null(external_csv) || !nzchar(external_csv) || !file.exists(external_csv)) {
    return(NULL)
  }
  External <- data.frame(data.table::fread(file = external_csv), check.names = FALSE)
  RN <- External[[1]]
  External <- External[, -1, drop = FALSE]
  row.names(External) <- RN
  n_ref <- ncol(ref_data)
  if (ncol(External) == n_ref - 1L) {
    colnames(External) <- colnames(ref_data)[seq_len(n_ref)[-1L]]
  } else if (ncol(External) == n_ref) {
    colnames(External) <- colnames(ref_data)
  } else {
    warning(sprintf(
      "External_Validation_Data.csv has %d columns after row-id drop; expected %d (no flag) or %d (full). Skipping external validation.",
      ncol(External), n_ref - 1L, n_ref
    ))
    return(NULL)
  }
  External$class <- factor(as.character(External$class), levels = levels(ref_data$class))
  External
}

find_external_validation_csv <- function(start_dirs, max_up = 8L) {
  hit <- env_csv_override(c("EXTERNAL_VALIDATION_CSV", "DEUTERATION_EXTERNAL_CSV"))
  if (!is.null(hit)) return(hit)

  search_from <- unique(normalizePath(start_dirs, winslash = "/", mustWork = FALSE))
  for (root in search_from) {
    cur <- root
    for (i in 0:max_up) {
      candidate <- file.path(cur, "Data", "External_Validation_Data.csv")
      if (file.exists(candidate)) {
        return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  NULL
}

training_csv <- find_training_csv(c(script_dir, getwd()))
if (is.null(training_csv)) {
  stop(
    "Could not find Data/Training_Data.csv. Set TRAINING_CSV to an explicit file path.",
    call. = FALSE
  )
}

message(sprintf("Loading data from: %s", training_csv))
ensure_kennard_stone_helpers(re_source = TRUE)
data <- load_training_data(training_csv)

external_csv <- find_external_validation_csv(c(script_dir, getwd()))
external_data <- load_external_validation_aligned(external_csv, data)
if (!is.null(external_data)) {
  message(sprintf("Loaded external validation data from: %s", external_csv))
} else {
  message("External validation data not found/aligned; external_accuracy will be NA.")
}

split_method <- "random"      # train/holdout: random | similarity | kennard_stone
initial_sampling <- "random"  # subset draw: random | stratified_min
min_per_class <- 5            # min rows per class when initial_sampling = stratified_min
balance_method <- "none"      # train balancing: none | smote | class_weight

output_dir <- Sys.getenv("DATA_SIZE_OUTPUT_DIR", unset = "")
if (!nzchar(output_dir)) {
  output_dir <- script_dir
}
output_dir <- ensure_dir_exists(output_dir, what = "output_dir")

for (nm in c(
  "write_progressive_csv", "save_progressive_csv_bundle",
  "init_inprogress_session", "load_inprogress_accumulators"
)) {
  if (exists(nm, envir = .GlobalEnv, inherits = FALSE)) {
    rm(list = nm, envir = .GlobalEnv)
  }
}

res <- evaluate_progressive_data_size(
  data = data,
  external_data = external_data,
  split_method = split_method,
  initial_sampling = initial_sampling,
  min_per_class = min_per_class,
  balance_method = balance_method,
  start_size = 20,                    # first subset size
  max_size = 50,                      # last subset size
  step_add = 5,                       # size step between evaluations
  train_p = 0.75,                     # holdout fraction (random split)
  n_repeats = 100,                    # repeats per size
  ordinal = TRUE,                     # polr vs multinom
  n_candidates = 200L,                # formulas searched on repeat 1
  top_k_per_repeat = 10L,             # top formulas ranked per repeat
  output_dir = output_dir,
  write_inprogress_csv = TRUE,        # snapshot CSVs after each size
  seed = 1
)

message("Progressive run completed. Summary metrics (best model per size, averaged over repeats):")
print(res$metrics_summary)
message("Best models selected by top-10 stability:")
print(res$best_models)
print(res$plot)

tryCatch(
  save_progressive_outputs(
    res = res,
    output_dir = output_dir,
    split_method = split_method,
    initial_sampling = initial_sampling,
    balance_method = balance_method
  ),
  error = function(e) {
    warning(
      "Final CSV export failed: ", conditionMessage(e),
      "\nRecover with: finalize_progressive_inprogress_to_final(output_dir, split_method, initial_sampling, balance_method)",
      call. = FALSE
    )
  }
)
