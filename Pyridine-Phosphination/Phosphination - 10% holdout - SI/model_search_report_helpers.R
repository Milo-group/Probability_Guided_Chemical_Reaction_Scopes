# PDF helpers for Model Search reports (used by Generate_Model_Search_Reports.R).
# Expects `report_dir` to be set in the calling environment before sourcing.

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush.console()
}

`%||%` <- function(x, y) if (is.null(x)) y else x

place_label <- function(i) {
  ord <- c("1st", "2nd", "3rd", paste0(4:10, "th"))[i]
  paste0(i, ". ", ord, " Place")
}

display_formula <- function(formula_chr) {
  gsub("`", "", as.character(formula_chr))
}

clean_var_name <- function(name) {
  gsub("`", "", as.character(name))
}

formula_predictor_names <- function(formula_chr) {
  vars <- all.vars(as.formula(formula_chr))
  vars[vars != "class"]
}

ensure_numeric_predictors <- function(data, formula_chr) {
  out <- data
  for (nm in formula_predictor_names(formula_chr)) {
    col <- clean_var_name(nm)
    if (col %in% names(out) && !is.numeric(out[[col]])) {
      out[[col]] <- as.numeric(as.character(out[[col]]))
    }
  }
  out
}

fit_candidate_model <- function(formula_chr, train_data, ordinal_flag) {
  data_fit <- ensure_numeric_predictors(train_data, formula_chr)
  if (ordinal_flag) {
    fit_polr(formula = formula_chr, data = data_fit)
  } else {
    nnet::multinom(formula_chr, data = data_fit, maxit = 2000, trace = FALSE)
  }
}

fit_models_with_accuracy <- function(models_df, train_data, test_data, ordinal_flag) {
  n <- nrow(models_df)
  train_acc <- numeric(n)
  test_acc <- numeric(n)

  for (i in seq_len(n)) {
    test.form <- models_df[i, 1][[1]]
    model <- tryCatch(
      fit_candidate_model(test.form, train_data, ordinal_flag),
      error = function(e) NULL
    )

    if (is.null(model)) {
      train_acc[i] <- NA_real_
      test_acc[i] <- NA_real_
    } else {
      mi_train <- mod.info(model, train_data, FALSE, TRUE)
      mi_test <- mod.info(model, test_data, FALSE, FALSE)
      train_acc[i] <- round(mi_train$accuracy, 2)
      test_acc[i] <- round(mi_test$accuracy, 2)
    }
  }

  cbind(
    models_df,
    `Training Accuracy` = train_acc,
    `Test Accuracy` = test_acc
  )
}

print_text_page <- function(title, lines, cex = 0.85) {
  plot.new()
  par(mar = c(1, 1, 1, 1))
  if (nzchar(title)) {
    text(0.5, 0.97, title, cex = 1.4, font = 2)
    y_top <- 0.9
  } else {
    y_top <- 0.97
  }
  text(0.05, y_top, paste(lines, collapse = "\n"), adj = c(0, 1), cex = cex, family = "mono")
}

format_toc_entry <- function(label, page, width = 90) {
  page_chr <- as.character(page)
  dots_n <- max(1L, width - nchar(label) - nchar(page_chr))
  paste0(label, paste(rep(".", dots_n), collapse = ""), page_chr)
}

format_plain_table <- function(df) {
  cols <- names(df)
  col_widths <- vapply(
    cols,
    function(col) max(nchar(col), max(nchar(as.character(df[[col]])), 0L)),
    integer(1)
  )
  fmt_row <- function(values) {
    paste(
      mapply(
        function(val, w) format(as.character(val), width = w, justify = "left"),
        values,
        col_widths
      ),
      collapse = "  "
    )
  }
  header <- fmt_row(cols)
  rows <- apply(df, 1, fmt_row)
  c(header, rows)
}

build_coef_table <- function(model, formula_chr, ordinal_flag) {
  pred_names <- formula_predictor_names(formula_chr)
  pred_clean <- vapply(pred_names, clean_var_name, character(1))
  expected_n <- length(pred_clean) + 1L

  if (ordinal_flag) {
    slopes <- coef(model)
    slope_vals <- vapply(pred_clean, function(pn) {
      idx <- match(pn, names(slopes))
      if (is.na(idx)) NA_real_ else as.numeric(slopes[idx])
    }, numeric(1))
    df <- data.frame(
      Term = c(pred_clean, "(Intercept)"),
      Estimate = c(slope_vals, as.numeric(model$zeta[1])),
      stringsAsFactors = FALSE
    )
  } else {
    cf <- coef(model)
    if (is.matrix(cf)) cf <- cf[1, , drop = TRUE]
    raw <- data.frame(
      Term = names(cf),
      Estimate = as.numeric(cf),
      stringsAsFactors = FALSE
    )
    keep <- raw$Term %in% c("(Intercept)", pred_clean)
    df <- raw[keep, , drop = FALSE]
    intercept_row <- df[df$Term == "(Intercept)", , drop = FALSE]
    pred_rows <- df[match(pred_clean, df$Term), , drop = FALSE]
    df <- rbind(intercept_row, pred_rows)
  }

  if (nrow(df) != expected_n) {
    warning(
      sprintf(
        "Coefficient count (%d) does not match formula predictors + intercept (%d): %s",
        nrow(df),
        expected_n,
        display_formula(formula_chr)
      ),
      call. = FALSE
    )
  }
  df
}

format_coef_lines <- function(model, formula_chr, ordinal_flag) {
  if (is.null(coef(model))) return("(coefficients unavailable)")
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Please install 'knitr' to format model coefficients in PDF reports.")
  }
  capture.output(knitr::kable(build_coef_table(model, formula_chr, ordinal_flag), row.names = FALSE))
}

format_model_stats_lines <- function(model, data, formula_chr, ordinal_flag) {
  mi <- mod.info(model, data, FALSE, FALSE)
  c(
    "Table:",
    "Full Model Stats - Overall Accuracy and Pseudo-R2",
    "",
    paste(sprintf("%-12s", "Accuracy"), "McFadden_R2"),
    paste(sprintf("%-12s", paste0(round(mi$accuracy, 2), "%")), mi$McFadden_R2),
    "",
    "Model Coefficients",
    format_coef_lines(model, formula_chr, ordinal_flag)
  )
}

write_model_search_report <- function(train_data,
                                      test_data,
                                      models,
                                      ordinal,
                                      report_title,
                                      output_file) {
  if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

  log_msg(sprintf("Writing Model Search report: %s", output_file))
  log_msg("  Computing train/test accuracies for top models ...")
  summary_table <- fit_models_with_accuracy(models, train_data, test_data, ordinal)
  mcf_col <- grep("McFadden", names(models), value = TRUE)[1]

  search_display <- data.frame(
    formula = vapply(models[[1]], display_formula, character(1)),
    `McFadden R2` = models[[mcf_col]],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  results_summary <- data.frame(
    Model = vapply(seq_len(nrow(summary_table)), place_label, character(1)),
    `Training Accuracy` = summary_table$`Training Accuracy`,
    `Test Accuracy` = summary_table$`Test Accuracy`,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  place_names <- vapply(
    seq_len(nrow(models)),
    function(i) gsub("^\\d+\\. ", "", place_label(i)),
    character(1)
  )

  render_pdf <- function(path, page_map = NULL, include_toc = FALSE) {
    if (is.null(page_map)) page_map <- list()
    page <- 0L

    next_page <- function(label = NULL) {
      page <<- page + 1L
      if (!is.null(label)) page_map[[label]] <<- page
      page
    }

    pdf(path, width = 11, height = 8.5)

    next_page("front")
    front_lines <- c(report_title, "", "Table of Contents")
    if (include_toc) {
      front_lines <- c(
        front_lines,
        format_toc_entry("Model Search", page_map[["model_search"]]),
        format_toc_entry("Section Results Summary", page_map[["results_summary"]]),
        vapply(
          place_names,
          function(nm) format_toc_entry(nm, page_map[[paste0("place_", nm)]]),
          character(1)
        ),
        "",
        "Model Search",
        format_plain_table(search_display),
        "",
        "Section Results Summary",
        format_plain_table(results_summary)
      )
    } else {
      page_map[["model_search"]] <- page
      page_map[["results_summary"]] <- page
      front_lines <- c(
        front_lines,
        "",
        "Model Search",
        format_plain_table(search_display),
        "",
        "Section Results Summary",
        format_plain_table(results_summary)
      )
    }
    print_text_page("", front_lines, cex = 0.72)

    fitted_models <- vector("list", nrow(models))
    for (i in seq_len(nrow(models))) {
      fitted_models[[i]] <- fit_candidate_model(models[i, 1][[1]], train_data, ordinal)
    }

    for (i in seq_len(nrow(models))) {
      conf_label <- place_label(i)
      place_name <- place_names[i]
      test.form <- models[i, 1][[1]]
      model <- fitted_models[[i]]

      next_page(paste0("place_", place_name))
      print_text_page(
        place_name,
        format_model_stats_lines(model, train_data, test.form, ordinal),
        cex = 0.8
      )

      model.info <- mod.info(model, train_data, FALSE, FALSE)
      next_page()
      print(
        ct_plot(
          model.info$class.table,
          plot.title = "Training Set",
          conformation = conf_label
        )$plot
      )
      next_page()
      print(
        prob.heatmap(
          model,
          train_data,
          plot.title = "Training Set",
          conformation = conf_label
        )
      )

      model.info <- mod.info(model, test_data, FALSE, FALSE)
      next_page()
      print(
        ct_plot(
          model.info$class.table,
          plot.title = "Test Set",
          conformation = conf_label
        )$plot
      )
      next_page()
      print(
        prob.heatmap(
          model,
          test_data,
          plot.title = "Test Set",
          conformation = conf_label
        )
      )
    }

    dev.off()
    page_map
  }

  log_msg("  Pass 1/2: measuring page layout ...")
  page_map <- render_pdf(tempfile(fileext = ".pdf"), include_toc = FALSE)

  for (i in seq_along(place_names)) {
    key <- paste0("place_", place_names[i])
    if (is.null(page_map[[key]])) {
      page_map[[key]] <- page_map[[paste0("place_", i)]] %||% (2L + (i - 1L) * 5L)
    }
  }

  log_msg("  Pass 2/2: rendering report with table of contents ...")
  temp_pdf <- tempfile(fileext = ".pdf")
  render_pdf(temp_pdf, page_map = page_map, include_toc = TRUE)

  if (!file.copy(temp_pdf, output_file, overwrite = TRUE)) {
    stop("Failed to write report PDF to: ", output_file)
  }
  unlink(temp_pdf)
  log_msg(sprintf("Model Search report written: %s", output_file))
}
