#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tools/benchmarks/benchmark_rcpp_cv.R"
script_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

trailing <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- trailing[startsWith(trailing, prefix)]
  if (length(hit)) sub(prefix, "", hit[[1]], fixed = TRUE) else default
}

base_ref <- arg_value("base-ref", "origin/dev")
first_pass_ref <- arg_value("first-pass-ref", "194373e")
out_dir <- normalizePath(arg_value("out-dir", script_dir), mustWork = FALSE)
reps <- as.integer(arg_value("reps", "3"))
quick <- "--quick" %in% trailing
threshold <- as.numeric(arg_value("threshold", "2"))
incremental_speed_threshold <- as.numeric(arg_value("incremental-speed-threshold", "1.10"))
rss_reduction_threshold <- as.numeric(arg_value("rss-reduction-threshold", "0.10"))
e2e_regression_threshold <- as.numeric(arg_value("e2e-regression-threshold", "0.95"))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!is.finite(reps) || reps < 1) stop("`--reps` must be a positive integer.")
if (!file.exists("/usr/bin/time")) stop("`/usr/bin/time` is required for RSS memory evidence.")

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
log_dir <- normalizePath(
  arg_value("log-dir", file.path(tempdir(), sprintf("rcpp_cv_time_logs_%s", timestamp))),
  mustWork = FALSE
)
results_csv <- file.path(out_dir, sprintf("rcpp_cv_benchmark_%s.csv", timestamp))
report_md <- file.path(out_dir, sprintf("rcpp_cv_benchmark_report_%s.md", timestamp))
latest_csv <- file.path(out_dir, "rcpp_cv_benchmark_latest.csv")
latest_report <- file.path(out_dir, "rcpp_cv_benchmark_report_latest.md")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

run_cmd <- function(command, args, stdout = TRUE, stderr = TRUE) {
  status <- system2(command, args, stdout = stdout, stderr = stderr)
  code <- attr(status, "status")
  if (is.null(code)) code <- 0L
  if (!identical(as.integer(code), 0L)) {
    if (is.character(status) && length(status)) {
      cat(paste(status, collapse = "\n"), "\n", file = stderr())
    }
    stop(sprintf("Command failed: %s %s", command, paste(args, collapse = " ")))
  }
  invisible(status)
}

copy_worktree_tracked <- function(repo, dest) {
  files <- system2("git", c("-C", repo, "ls-files"), stdout = TRUE)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  for (rel in files) {
    src <- file.path(repo, rel)
    dst <- file.path(dest, rel)
    if (!file.exists(src)) next
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, dst, overwrite = TRUE, copy.date = TRUE)
  }
}

archive_ref <- function(repo, ref, dest) {
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  tar_file <- tempfile(fileext = ".tar")
  on.exit(unlink(tar_file), add = TRUE)
  run_cmd("git", c("-C", repo, "archive", "--format=tar", paste0("--output=", tar_file), ref))
  untar(tar_file, exdir = dest)
}

install_pkg <- function(src, lib, label) {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_dir, sprintf("install_%s_%s.log", label, timestamp))
  status <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-multiarch", "--with-keep.source", paste0("--library=", lib), src),
    stdout = log_file,
    stderr = log_file
  )
  code <- attr(status, "status")
  if (is.null(code)) code <- 0L
  if (!identical(as.integer(code), 0L)) {
    stop(sprintf("Package install failed for %s; see %s", label, log_file))
  }
}

scenario_id <- function(estimator, link, size, pop_totals, est_method) {
  sprintf(
    "%s|%s|%s|pop=%s|est=%s",
    estimator,
    link,
    size,
    pop_totals,
    if (is.na(est_method)) "none" else est_method
  )
}

make_scenarios <- function() {
  links <- c("logit", "probit", "cloglog")
  rows <- list()
  k <- 0L
  for (link in links) {
    for (size_name in c("moderate", "high")) {
      for (use_totals in c(FALSE, TRUE)) {
        k <- k + 1L
        rows[[k]] <- data.frame(
          estimator = "direct_cv",
          link = link,
          size = size_name,
          pop_totals = use_totals,
          est_method = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
    for (est_method in c("mle", "gee")) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        estimator = "end_to_end",
        link = link,
        size = "moderate",
        pop_totals = FALSE,
        est_method = est_method,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$scenario_id <- mapply(
    scenario_id,
    out$estimator,
    out$link,
    out$size,
    out$pop_totals,
    out$est_method,
    USE.NAMES = FALSE
  )
  out
}

bench_worker <- tempfile(fileext = ".R")
writeLines(c(
  "args <- commandArgs(trailingOnly = TRUE)",
  "lib <- args[[1]]",
  "out <- args[[2]]",
  "reps <- as.integer(args[[3]])",
  "quick <- identical(args[[4]], 'TRUE')",
  "estimator <- args[[5]]",
  "link <- args[[6]]",
  "size_name <- args[[7]]",
  "use_totals <- identical(args[[8]], 'TRUE')",
  "est_method <- if (identical(args[[9]], 'none')) NA_character_ else args[[9]]",
  ".libPaths(c(lib, .libPaths()))",
  "suppressPackageStartupMessages(library(nonprobsvy))",
  "suppressPackageStartupMessages(library(survey))",
  "",
  "time_expr <- function(expr, reps) {",
  "  expr <- substitute(expr)",
  "  env <- parent.frame()",
  "  timings <- numeric(reps)",
  "  value <- NULL",
  "  for (i in seq_len(reps)) {",
  "    gc(FALSE)",
  "    timings[[i]] <- system.time(value <- eval(expr, envir = env))[['elapsed']]",
  "  }",
  "  list(time = timings, value = value)",
  "}",
  "",
  "make_direct_data <- function(seed, n, p) {",
  "  set.seed(seed)",
  "  X <- scale(matrix(rnorm(n * p), nrow = n, ncol = p))",
  "  beta <- c(0.6, -0.45, 0.3, rep(0, p - 3))",
  "  eta <- drop(X %*% beta)",
  "  R <- rbinom(n, size = 1, prob = plogis(eta))",
  "  list(X = X, R = R, weights = runif(n, 0.8, 1.4), pop_totals = colSums(X) + seq_len(p) / 100)",
  "}",
  "",
  "run_direct <- function(link, size_name, use_totals, reps) {",
  "  grid <- if (size_name == 'high') c(n = if (quick) 220 else 420, p = if (quick) 18 else 36) else c(n = if (quick) 160 else 260, p = if (quick) 10 else 18)",
  "  dat <- make_direct_data(1000 + grid[['n']] + grid[['p']], grid[['n']], grid[['p']])",
  "  pop <- if (use_totals) dat$pop_totals else NULL",
  "  timed <- time_expr({",
  "    set.seed(3001)",
  "    nonprobsvy:::cv_nonprobsvy_rcpp(",
  "      X = dat$X, R = dat$R, weights_X = dat$weights, method_selection = link,",
  "      gee_h_fun = if (use_totals) 1 else 2, maxit = if (quick) 45 else 80, eps = 1e-4,",
  "      lambda_min = 0.01, nlambda = if (quick) 4 else 8, nfolds = if (quick) 2 else 3,",
  "      penalty = 'lasso', a = 3, pop_totals = pop, verbose = FALSE)",
  "  }, reps)",
  "  val <- timed$value",
  "  list(times = timed$time, summary = list(",
  "    theta = as.numeric(val$theta_est), selected = as.integer(val$theta_selected),",
  "    lambda = as.numeric(val$lambda), cv_error = as.numeric(val$cv_error)))",
  "}",
  "",
  "make_e2e_data <- function(seed, n_a, n_b, p) {",
  "  set.seed(seed)",
  "  x_a <- matrix(rnorm(n_a * p), n_a, p)",
  "  x_b <- matrix(rnorm(n_b * p), n_b, p)",
  "  beta <- c(0.5, -0.35, 0.25, rep(0, p - 3))",
  "  y_a <- 1 + drop(x_a %*% beta) + rnorm(n_a, sd = 0.5)",
  "  y_b <- 1 + drop(x_b %*% beta) + rnorm(n_b, sd = 0.5)",
  "  nons <- as.data.frame(x_a)",
  "  prob <- as.data.frame(x_b)",
  "  names(nons) <- names(prob) <- paste0('x', seq_len(p))",
  "  nons$y <- y_a",
  "  prob$y <- y_b",
  "  prob$w <- runif(n_b, 0.8, 1.6)",
  "  list(nons = nons, prob = prob)",
  "}",
  "",
  "run_e2e <- function(link, est_method, use_totals, reps) {",
  "  p <- if (quick) 6 else 8",
  "  dat <- make_e2e_data(5000 + p + nchar(link) + nchar(est_method) + 100 * use_totals, if (quick) 100 else 160, if (quick) 70 else 110, p)",
  "  selection <- reformulate(paste0('x', seq_len(p)))",
  "  control_selection <- control_sel(est_method = est_method, gee_h_fun = 1, maxit = if (quick) 45 else 70, nfolds = 2, nlambda = if (quick) 4 else 6, penalty = 'lasso')",
  "  timed <- time_expr({",
  "    set.seed(7001)",
  "    if (use_totals) {",
  "      mm <- model.matrix(selection, dat$prob)",
  "      pop_totals <- c(colSums(mm * dat$prob$w), y = sum(dat$prob$y * dat$prob$w))",
  "      suppressWarnings(nonprob(",
  "        selection = selection, target = ~y, pop_totals = pop_totals, data = dat$nons,",
  "        method_selection = link, control_inference = control_inf(vars_selection = TRUE),",
  "        control_selection = control_selection, se = FALSE))",
  "    } else {",
  "      svy <- svydesign(ids = ~1, weights = ~w, data = dat$prob)",
  "      suppressWarnings(nonprob(",
  "        selection = selection, target = ~y, svydesign = svy, data = dat$nons,",
  "        method_selection = link, control_inference = control_inf(vars_selection = TRUE),",
  "        control_selection = control_selection, se = FALSE))",
  "    }",
  "  }, reps)",
  "  val <- timed$value",
  "  list(times = timed$time, summary = list(",
  "    mean = as.numeric(val$output$mean), selected = names(coef(val)$coef_sel[, 1])))",
  "}",
  "",
  "result <- if (identical(estimator, 'direct_cv')) {",
  "  run_direct(link, size_name, use_totals, reps)",
  "} else {",
  "  run_e2e(link, est_method, use_totals, reps)",
  "}",
  "saveRDS(list(",
  "  estimator = estimator, link = link, size = size_name, pop_totals = use_totals,",
  "  est_method = est_method, times = result$times, summary = result$summary), out)"
), bench_worker)

prepare_source <- function(label, ref = NULL) {
  src <- file.path(tempdir(), paste0("nonprobsvy-src-", label, "-", timestamp))
  if (is.null(ref)) {
    copy_worktree_tracked(repo_root, src)
  } else {
    archive_ref(repo_root, ref, src)
  }
  src
}

parse_max_rss <- function(log_file) {
  lines <- readLines(log_file, warn = FALSE)
  hit <- grep("maximum resident set size", lines, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(sub("^\\s*([0-9]+).*$", "\\1", hit[[length(hit)]]))
}

run_worker <- function(label, lib, scenario, index) {
  out <- tempfile(fileext = ".rds")
  log_file <- file.path(log_dir, sprintf("time_%s_%02d_%s.log", label, index, timestamp))
  est_method_arg <- if (is.na(scenario$est_method)) "none" else scenario$est_method
  worker_args <- c(
    "-l",
    file.path(R.home("bin"), "Rscript"),
    bench_worker,
    lib,
    out,
    as.character(reps),
    as.character(quick),
    scenario$estimator,
    scenario$link,
    scenario$size,
    as.character(scenario$pop_totals),
    est_method_arg
  )

  status <- system2("/usr/bin/time", worker_args, stdout = log_file, stderr = log_file)
  code <- attr(status, "status")
  if (is.null(code)) code <- 0L
  if (!identical(as.integer(code), 0L)) {
    stop(sprintf("Benchmark worker failed for %s / %s; see %s", label, scenario$scenario_id, log_file))
  }

  row <- readRDS(out)
  row$scenario_id <- scenario$scenario_id
  row$max_rss_bytes <- parse_max_rss(log_file)
  row
}

run_impl <- function(label, src, scenarios) {
  lib <- file.path(tempdir(), paste0("nonprobsvy-lib-", label, "-", timestamp))
  install_pkg(src, lib, label)
  lapply(seq_len(nrow(scenarios)), function(i) {
    message(sprintf("Running %s scenario %d/%d: %s", label, i, nrow(scenarios), scenarios$scenario_id[[i]]))
    run_worker(label, lib, scenarios[i, ], i)
  })
}

summarize_one <- function(label, rows) {
  do.call(rbind, lapply(rows, function(row) {
    data.frame(
      implementation = label,
      scenario_id = row$scenario_id,
      estimator = row$estimator,
      link = row$link,
      size = row$size,
      pop_totals = row$pop_totals,
      est_method = if (is.na(row$est_method)) NA_character_ else row$est_method,
      median_sec = median(row$times),
      min_sec = min(row$times),
      max_sec = max(row$times),
      max_rss_bytes = row$max_rss_bytes,
      stringsAsFactors = FALSE
    )
  }))
}

numeric_summary <- function(reference, current) {
  if (!identical(reference$estimator, current$estimator)) stop("Scenario mismatch.")
  if (identical(reference$estimator, "direct_cv")) {
    theta_diff <- max(abs(reference$summary$theta - current$summary$theta))
    cv_diff <- max(abs(reference$summary$cv_error - current$summary$cv_error))
    lambda_diff <- abs(reference$summary$lambda - current$summary$lambda)
    selected_same <- identical(reference$summary$selected, current$summary$selected)
    list(max_abs_diff = max(theta_diff, lambda_diff), cv_error_max_abs_diff = cv_diff, selected_identical = selected_same)
  } else {
    mean_diff <- max(abs(reference$summary$mean - current$summary$mean))
    selected_same <- identical(reference$summary$selected, current$summary$selected)
    list(max_abs_diff = mean_diff, cv_error_max_abs_diff = NA_real_, selected_identical = selected_same)
  }
}

add_metrics <- function(wide, summary, suffix) {
  idx <- match(wide$scenario_id, summary$scenario_id)
  for (metric in c("median_sec", "min_sec", "max_sec", "max_rss_bytes")) {
    wide[[paste0(metric, "_", suffix)]] <- summary[[metric]][idx]
  }
  wide
}

add_numeric_metrics <- function(wide, reference_rows, current_rows, suffix) {
  idx_ref <- match(wide$scenario_id, vapply(reference_rows, `[[`, character(1), "scenario_id"))
  idx_cur <- match(wide$scenario_id, vapply(current_rows, `[[`, character(1), "scenario_id"))
  num <- Map(numeric_summary, reference_rows[idx_ref], current_rows[idx_cur])
  wide[[paste0("max_abs_diff_vs_", suffix)]] <- vapply(num, `[[`, numeric(1), "max_abs_diff")
  wide[[paste0("cv_error_max_abs_diff_vs_", suffix)]] <- vapply(num, `[[`, numeric(1), "cv_error_max_abs_diff")
  wide[[paste0("selected_identical_vs_", suffix)]] <- vapply(num, `[[`, logical(1), "selected_identical")
  wide
}

fmt_sec <- function(x) sprintf("%.4fs", x)
fmt_x <- function(x) sprintf("%.3fx", x)
fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt_mb <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f MB", x / 1024^2))
fmt_diff <- function(x) format(x, scientific = TRUE, digits = 4)

scenarios <- make_scenarios()

message("Preparing and installing baseline: ", base_ref)
base_src <- prepare_source("baseline_dev", base_ref)
base_rows <- run_impl("baseline_dev", base_src, scenarios)

message("Preparing and installing first pass: ", first_pass_ref)
first_src <- prepare_source("first_pass", first_pass_ref)
first_rows <- run_impl("first_pass", first_src, scenarios)

message("Preparing and installing current worktree")
current_src <- prepare_source("current", NULL)
current_rows <- run_impl("current", current_src, scenarios)

base_summary <- summarize_one("baseline_dev", base_rows)
first_summary <- summarize_one("first_pass", first_rows)
current_summary <- summarize_one("current", current_rows)

wide <- scenarios[, c("scenario_id", "estimator", "link", "size", "pop_totals", "est_method")]
wide <- add_metrics(wide, base_summary, "baseline")
wide <- add_metrics(wide, first_summary, "first_pass")
wide <- add_metrics(wide, current_summary, "current")
wide <- add_numeric_metrics(wide, base_rows, current_rows, "baseline")
wide <- add_numeric_metrics(wide, first_rows, current_rows, "first_pass")

wide$speedup_vs_baseline <- wide$median_sec_baseline / wide$median_sec_current
wide$speedup_vs_first_pass <- wide$median_sec_first_pass / wide$median_sec_current
wide$rss_ratio_vs_baseline <- wide$max_rss_bytes_baseline / wide$max_rss_bytes_current
wide$rss_ratio_vs_first_pass <- wide$max_rss_bytes_first_pass / wide$max_rss_bytes_current
wide$rss_reduction_vs_first_pass <- 1 - wide$max_rss_bytes_current / wide$max_rss_bytes_first_pass

direct <- wide[wide$estimator == "direct_cv", , drop = FALSE]
e2e <- wide[wide$estimator == "end_to_end", , drop = FALSE]
median_direct_speedup_vs_baseline <- median(direct$speedup_vs_baseline)
median_direct_speedup_vs_first <- median(direct$speedup_vs_first_pass)
median_direct_rss_reduction_vs_first <- median(direct$rss_reduction_vs_first_pass, na.rm = TRUE)
numeric_pass_first <- all(wide$max_abs_diff_vs_first_pass <= 1e-6) && all(wide$selected_identical_vs_first_pass)
numeric_pass_baseline <- all(wide$max_abs_diff_vs_baseline <= 1e-6) && all(wide$selected_identical_vs_baseline)
end_to_end_regression_free_first <- all(e2e$speedup_vs_first_pass >= e2e_regression_threshold)
incremental_runtime_pass <- median_direct_speedup_vs_first >= incremental_speed_threshold
incremental_memory_pass <- is.finite(median_direct_rss_reduction_vs_first) &&
  median_direct_rss_reduction_vs_first >= rss_reduction_threshold
incremental_gate_pass <- numeric_pass_first &&
  end_to_end_regression_free_first &&
  (incremental_runtime_pass || incremental_memory_pass)

write.csv(wide, results_csv, row.names = FALSE)
file.copy(results_csv, latest_csv, overwrite = TRUE)

scenario_lines <- sprintf(
  paste0(
    "- `%s`: time dev `%s`, first `%s`, current `%s`; ",
    "speedup current/dev `%s`, current/first `%s`; ",
    "RSS dev `%s`, first `%s`, current `%s`; ",
    "RSS reduction current/first `%s`; ",
    "max fitted diff vs first `%s`, CV-error diff vs first `%s`, selected vs first `%s`"
  ),
  wide$scenario_id,
  fmt_sec(wide$median_sec_baseline),
  fmt_sec(wide$median_sec_first_pass),
  fmt_sec(wide$median_sec_current),
  fmt_x(wide$speedup_vs_baseline),
  fmt_x(wide$speedup_vs_first_pass),
  fmt_mb(wide$max_rss_bytes_baseline),
  fmt_mb(wide$max_rss_bytes_first_pass),
  fmt_mb(wide$max_rss_bytes_current),
  fmt_pct(wide$rss_reduction_vs_first_pass),
  fmt_diff(wide$max_abs_diff_vs_first_pass),
  ifelse(is.na(wide$cv_error_max_abs_diff_vs_first_pass), "NA", fmt_diff(wide$cv_error_max_abs_diff_vs_first_pass)),
  wide$selected_identical_vs_first_pass
)

report <- c(
  "# Rcpp CV Benchmark Report",
  "",
  sprintf("- Timestamp: `%s`", timestamp),
  sprintf("- Baseline ref: `%s`", base_ref),
  sprintf("- First-pass ref: `%s`", first_pass_ref),
  sprintf("- Current ref: worktree at `%s`", repo_root),
  sprintf("- Repetitions per scenario: `%s`", reps),
  sprintf("- Quick mode: `%s`", quick),
  sprintf("- Required median direct-CV speedup versus dev: `%.2fx`", threshold),
  sprintf("- Incremental speed threshold versus first pass: `%.2fx`", incremental_speed_threshold),
  sprintf("- Incremental max-RSS reduction threshold versus first pass: `%s`", fmt_pct(rss_reduction_threshold)),
  "",
  "## Gate",
  "",
  sprintf("- Median direct-CV speedup versus dev: `%.3fx`", median_direct_speedup_vs_baseline),
  sprintf("- Median direct-CV speedup versus first pass: `%.3fx`", median_direct_speedup_vs_first),
  sprintf("- Median direct-CV max-RSS reduction versus first pass: `%s`", fmt_pct(median_direct_rss_reduction_vs_first)),
  sprintf("- Numerical equivalence versus first pass (`max_abs_diff <= 1e-6` and identical selected sets): `%s`", numeric_pass_first),
  sprintf("- Numerical equivalence versus dev (`max_abs_diff <= 1e-6` and identical selected sets): `%s`", numeric_pass_baseline),
  sprintf("- End-to-end regression-free versus first pass (`speedup >= %.2fx`): `%s`", e2e_regression_threshold, end_to_end_regression_free_first),
  sprintf("- Incremental runtime pass: `%s`", incremental_runtime_pass),
  sprintf("- Incremental memory pass: `%s`", incremental_memory_pass),
  sprintf("- Overall incremental gate: `%s`", if (incremental_gate_pass) "PASS" else "FAIL"),
  "",
  "## Scenario Results",
  "",
  scenario_lines,
  "",
  "## Artifacts",
  "",
  sprintf("- CSV: `%s`", basename(results_csv)),
  sprintf("- Latest CSV copy: `%s`", basename(latest_csv))
)

writeLines(report, report_md)
file.copy(report_md, latest_report, overwrite = TRUE)

message("Wrote: ", results_csv)
message("Wrote: ", report_md)
message("Overall incremental gate: ", if (incremental_gate_pass) "PASS" else "FAIL")

if (!incremental_gate_pass) quit(status = 2L)
