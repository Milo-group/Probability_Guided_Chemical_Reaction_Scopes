# Kennard–Stone (KS) subset selection in feature space (class-wise wrapper).
# Used by Model Construction / Cross validation scripts and progressive data-size evaluation.

prepare_deuteration_training_data <- function(data_or_path) {
  if (is.character(data_or_path) && length(data_or_path) == 1L && file.exists(data_or_path)) {
    data <- data.frame(data.table::fread(data_or_path, check.names = FALSE))
  } else if (is.data.frame(data_or_path)) {
    data <- data_or_path
  } else {
    stop("provide a data.frame or path to Training_Data.csv", call. = FALSE)
  }
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

source_kennard_stone_helpers <- function(ks_path, local = FALSE) {
  ks_path <- normalizePath(ks_path, winslash = "/", mustWork = FALSE)
  if (!file.exists(ks_path)) {
    stop("Kennard–Stone helpers not found: ", ks_path, call. = FALSE)
  }
  source(ks_path, local = local)
  invisible(ks_path)
}

coerce_predictor_matrix <- function(data, feature_cols) {
  feature_cols <- intersect(as.character(feature_cols), names(data))
  feature_cols <- feature_cols[nzchar(feature_cols)]
  if (!length(feature_cols)) {
    return(matrix(numeric(0), nrow = nrow(data), ncol = 0L))
  }

  n <- nrow(data)
  cols <- lapply(feature_cols, function(nm) {
    suppressWarnings(as.numeric(data[[nm]]))
  })
  mat <- do.call(cbind, cols)
  if (is.null(mat)) {
    return(matrix(numeric(0), nrow = n, ncol = 0L))
  }
  storage.mode(mat) <- "double"
  colnames(mat) <- feature_cols

  bad <- !is.finite(mat)
  if (any(bad)) {
    for (j in seq_len(ncol(mat))) {
      col_j <- mat[, j]
      bad_j <- !is.finite(col_j)
      if (any(bad_j)) {
        fill <- stats::median(col_j[!bad_j], na.rm = TRUE)
        if (!is.finite(fill)) {
          fill <- 0
        }
        col_j[bad_j] <- fill
        mat[, j] <- col_j
      }
    }
  }
  mat
}

# Predictors used for Kennard–Stone distances (exclude response and row index).
ks_feature_cols <- function(data, class_col = "class") {
  cols <- setdiff(names(data), c(class_col, "flag"))
  cols <- intersect(cols, names(data))
  cols[nzchar(cols)]
}

kennard_stone <- function(X, n_select) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)

  if (n_select >= n) {
    return(seq_len(n))
  }

  D <- as.matrix(stats::dist(X, method = "euclidean"))

  max_pair <- which(D == max(D), arr.ind = TRUE)[1, ]
  selected <- unique(as.integer(max_pair))

  while (length(selected) < n_select) {
    remaining <- setdiff(seq_len(n), selected)
    min_dist <- apply(D[remaining, selected, drop = FALSE], 1, min)
    next_sample <- remaining[which.max(min_dist)]
    selected <- c(selected, next_sample)
  }

  selected
}

classwise_kennard_stone <- function(data, feature_cols, class_col, n_per_class) {
  selected_rows <- integer(0)

  for (cl in unique(data[[class_col]])) {
    idx <- which(data[[class_col]] == cl)
    n_select <- min(as.integer(n_per_class), length(idx))
    if (n_select <= 0L) {
      next
    }
    X_class <- coerce_predictor_matrix(data[idx, , drop = FALSE], feature_cols)
    selected_local <- kennard_stone(X_class, n_select)
    selected_rows <- c(selected_rows, idx[selected_local])
  }

  unique(selected_rows)
}

scale_predictors <- function(data, feature_cols, scale_center = NULL, scale_scale = NULL) {
  mat <- coerce_predictor_matrix(data, feature_cols)
  feature_cols <- colnames(mat)
  if (!length(feature_cols) || ncol(mat) == 0L || nrow(mat) == 0L) {
    stop("No numeric predictor columns available for scaling.", call. = FALSE)
  }

  reuse_scaling <- !is.null(scale_center) &&
    !is.null(scale_scale) &&
    length(scale_center) == ncol(mat) &&
    length(scale_scale) == ncol(mat)

  scaled <- if (reuse_scaling) {
    base::scale(mat, center = scale_center, scale = scale_scale)
  } else {
    base::scale(mat)
  }

  out <- data
  out[feature_cols] <- as.data.frame(scaled, stringsAsFactors = FALSE)
  attr(out, "ks_scale_center") <- attr(scaled, "scaled:center")
  attr(out, "ks_scale_scale") <- attr(scaled, "scaled:scale")
  attr(out, "ks_feature_cols") <- feature_cols
  out
}

# Train indices: up to n_per_class diverse representatives per class (features scaled globally).
build_kennard_stone_train_indices <- function(data,
                                              class_col = "class",
                                              feature_cols = NULL,
                                              n_per_class = NULL,
                                              train_fraction = NULL) {
  if (is.null(feature_cols)) {
    feature_cols <- ks_feature_cols(data, class_col)
  }
  if (!length(feature_cols)) {
    stop("No predictor columns found for Kennard–Stone sampling.", call. = FALSE)
  }
  if (!class_col %in% names(data)) {
    stop(sprintf("class column '%s' not found.", class_col), call. = FALSE)
  }

  tab <- table(data[[class_col]])
  if (length(tab) == 0L) {
    return(integer(0))
  }

  if (is.null(n_per_class)) {
    if (!is.null(train_fraction)) {
      n_per_class <- max(1L, floor(min(as.integer(tab)) * train_fraction))
    } else {
      n_per_class <- as.integer(min(tab))
    }
  }
  n_per_class <- max(1L, as.integer(n_per_class))

  df_scaled <- scale_predictors(data, feature_cols)
  feature_cols <- attr(df_scaled, "ks_feature_cols")
  if (is.null(feature_cols) || !length(feature_cols)) {
    feature_cols <- ks_feature_cols(data, class_col)
  }
  classwise_kennard_stone(
    data = df_scaled,
    feature_cols = feature_cols,
    class_col = class_col,
    n_per_class = n_per_class
  )
}

# Split a data frame into train / holdout using class-wise Kennard–Stone (cf. simi.sampler splits).
split_with_kennard_stone_sampler <- function(data,
                                             seed = NULL,
                                             target_train_n = NULL,
                                             n_per_class = NULL,
                                             min_per_class = 5L,
                                             class_col = "class") {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  stratified_fallback <- function(df) {
    cls <- sort(unique(as.character(df[[class_col]])))
    picked <- integer(0)
    for (cl in cls) {
      idx <- which(as.character(df[[class_col]]) == cl)
      if (length(idx) == 0) {
        next
      }
      take_n <- max(1L, floor(length(idx) * 0.7))
      take_n <- min(take_n, length(idx))
      picked <- c(picked, sample(idx, size = take_n))
    }
    picked <- unique(picked)
    if (length(picked) == 0L) {
      picked <- sample(seq_len(nrow(df)), size = max(1L, floor(0.7 * nrow(df))))
    }
    list(
      train_data = df[picked, , drop = FALSE],
      test_data = df[-picked, , drop = FALSE]
    )
  }

  n <- nrow(data)
  if (n < 2L) {
    return(stratified_fallback(data))
  }

  feature_cols <- ks_feature_cols(data, class_col)
  if (!length(feature_cols)) {
    warning("No feature columns for Kennard–Stone; using stratified fallback split.")
    return(stratified_fallback(data))
  }

  if (is.null(target_train_n)) {
    target_train_n <- max(1L, floor(0.75 * n))
  }
  target_train_n <- max(1L, min(as.integer(target_train_n), n - 1L))

  if (is.null(n_per_class)) {
    n_classes <- length(unique(data[[class_col]]))
    n_per_class <- max(as.integer(min_per_class), ceiling(target_train_n / n_classes))
  }

  selected_idx <- tryCatch(
    build_kennard_stone_train_indices(
      data = data,
      class_col = class_col,
      feature_cols = feature_cols,
      n_per_class = n_per_class
    ),
    error = function(e) {
      message(sprintf(
        "Kennard–Stone split failed (n=%d): %s. Using stratified fallback.",
        n, conditionMessage(e)
      ))
      NULL
    }
  )

  if (is.null(selected_idx) || length(selected_idx) == 0L) {
    return(stratified_fallback(data))
  }

  if (length(selected_idx) > target_train_n) {
    df_scaled <- scale_predictors(data[selected_idx, , drop = FALSE], feature_cols)
    local_keep <- kennard_stone(
      as.matrix(df_scaled[, feature_cols, drop = FALSE]),
      target_train_n
    )
    selected_idx <- selected_idx[local_keep]
  } else if (length(selected_idx) < target_train_n) {
    remaining <- setdiff(seq_len(n), selected_idx)
    need <- target_train_n - length(selected_idx)
    if (length(remaining) >= need) {
      pool <- data[c(selected_idx, remaining), , drop = FALSE]
      df_scaled <- scale_predictors(pool, feature_cols)
      rel_remaining <- seq_along(remaining) + length(selected_idx)
      X_rem <- as.matrix(df_scaled[rel_remaining, feature_cols, drop = FALSE])
      if (length(selected_idx) > 0L) {
        X_sel <- as.matrix(df_scaled[seq_along(selected_idx), feature_cols, drop = FALSE])
        D <- as.matrix(stats::dist(rbind(X_sel, X_rem), method = "euclidean"))
        min_to_sel <- apply(D[(length(selected_idx) + 1L):nrow(D), seq_len(length(selected_idx)), drop = FALSE], 1, min)
        ord <- order(min_to_sel, decreasing = TRUE, na.last = NA)
        extra <- remaining[ord[seq_len(min(need, length(ord)))]]
      } else {
        extra_local <- kennard_stone(X_rem, need)
        extra <- remaining[extra_local]
      }
      selected_idx <- c(selected_idx, extra)
    }
  }

  selected_idx <- unique(selected_idx)
  selected_idx <- selected_idx[selected_idx >= 1L & selected_idx <= n]
  if (length(selected_idx) == 0L || length(selected_idx) >= n) {
    return(stratified_fallback(data))
  }

  list(
    train_data = data[selected_idx, , drop = FALSE],
    test_data = data[-selected_idx, , drop = FALSE]
  )
}

# Progressive subset: class-wise KS (min_per_class each), then KS fill to subset_size.
sample_kennard_stone_indices <- function(data,
                                         subset_size,
                                         class_col = "class",
                                         min_per_class = 5L,
                                         feature_cols = NULL) {
  n_total <- nrow(data)
  subset_size <- as.integer(subset_size)[1L]
  if (is.na(subset_size) || subset_size <= 0L || n_total == 0L) {
    return(integer(0))
  }
  subset_size <- min(subset_size, n_total)

  if (is.null(feature_cols)) {
    feature_cols <- ks_feature_cols(data, class_col)
  }
  if (!length(feature_cols)) {
    stop("kennard_stone sampling requires predictor columns.", call. = FALSE)
  }

  df_scaled <- scale_predictors(data, feature_cols)
  feature_cols <- attr(df_scaled, "ks_feature_cols")
  if (is.null(feature_cols) || !length(feature_cols)) {
    feature_cols <- ks_feature_cols(data, class_col)
  }

  selected <- classwise_kennard_stone(
    data = df_scaled,
    feature_cols = feature_cols,
    class_col = class_col,
    n_per_class = min_per_class
  )

  selected <- unique(selected)
  if (length(selected) > subset_size) {
    sub_mat <- coerce_predictor_matrix(df_scaled[selected, , drop = FALSE], feature_cols)
    local_keep <- kennard_stone(sub_mat, subset_size)
    selected <- selected[local_keep]
  } else if (length(selected) < subset_size) {
    remaining <- setdiff(seq_len(n_total), selected)
    need <- subset_size - length(selected)
    if (length(remaining) > 0L && need > 0L) {
      take_n <- min(need, length(remaining))
      rem_mat <- coerce_predictor_matrix(df_scaled[remaining, , drop = FALSE], feature_cols)
      if (length(selected) > 0L) {
        sel_mat <- coerce_predictor_matrix(df_scaled[selected, , drop = FALSE], feature_cols)
        combined <- rbind(sel_mat, rem_mat)
        D <- as.matrix(stats::dist(combined, method = "euclidean"))
        n_sel <- nrow(sel_mat)
        min_to_sel <- apply(D[(n_sel + 1L):nrow(D), seq_len(n_sel), drop = FALSE], 1, min)
        ord <- order(min_to_sel, decreasing = TRUE, na.last = NA)
        extra <- remaining[ord[seq_len(take_n)]]
      } else {
        extra <- remaining[kennard_stone(rem_mat, take_n)]
      }
      selected <- c(selected, extra)
    }
  }

  sort(unique(selected))
}
