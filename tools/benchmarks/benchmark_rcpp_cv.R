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
out_dir <- normalizePath(arg_value("out-dir", script_dir), mustWork = FALSE)
reps <- as.integer(arg_value("reps", "3"))
quick <- "--quick" %in% trailing
threshold <- as.numeric(arg_value("threshold", "2"))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!is.finite(reps) || reps < 1) stop("`--reps` must be a positive integer.")

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
results_csv <- file.path(out_dir, sprintf("rcpp_cv_benchmark_%s.csv", timestamp))
report_md <- file.path(out_dir, sprintf("rcpp_cv_benchmark_report_%s.md", timestamp))
latest_csv <- file.path(out_dir, "rcpp_cv_benchmark_latest.csv")
latest_report <- file.path(out_dir, "rcpp_cv_benchmark_report_latest.md")

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
  log_file <- file.path(out_dir, sprintf("install_%s_%s.log", label, timestamp))
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

bench_worker <- tempfile(fileext = ".R")
writeLines(c(
  "args <- commandArgs(trailingOnly = TRUE)",
  "lib <- args[[1]]",
  "out <- args[[2]]",
  "reps <- as.integer(args[[3]])",
  "quick <- identical(args[[4]], 'TRUE')",
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
  "links <- c('logit', 'probit', 'cloglog')",
  "rows <- list()",
  "k <- 0L",
  "for (link in links) {",
  "  for (size_name in c('moderate', 'high')) {",
  "    for (use_totals in c(FALSE, TRUE)) {",
  "      k <- k + 1L",
  "      rows[[k]] <- c(list(estimator = 'direct_cv', link = link, size = size_name, pop_totals = use_totals, est_method = NA_character_), run_direct(link, size_name, use_totals, reps))",
  "    }",
  "  }",
  "  for (est_method in c('mle', 'gee')) {",
  "    k <- k + 1L",
  "    rows[[k]] <- c(list(estimator = 'end_to_end', link = link, size = 'moderate', pop_totals = FALSE, est_method = est_method), run_e2e(link, est_method, FALSE, reps))",
  "  }",
  "}",
  "saveRDS(rows, out)"
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

summarize_one <- function(label, rows) {
  do.call(rbind, lapply(seq_along(rows), function(i) {
    row <- rows[[i]]
    data.frame(
      implementation = label,
      scenario_id = sprintf(
        "%s|%s|%s|pop=%s|est=%s",
        row$estimator,
        row$link,
        row$size,
        row$pop_totals,
        if (is.na(row$est_method)) "none" else row$est_method
      ),
      estimator = row$estimator,
      link = row$link,
      size = row$size,
      pop_totals = row$pop_totals,
      est_method = if (is.na(row$est_method)) NA_character_ else row$est_method,
      median_sec = median(row$times),
      min_sec = min(row$times),
      max_sec = max(row$times),
      stringsAsFactors = FALSE
    )
  }))
}

numeric_summary <- function(base, current) {
  if (!identical(base$estimator, current$estimator)) stop("Scenario mismatch.")
  if (identical(base$estimator, "direct_cv")) {
    theta_diff <- max(abs(base$summary$theta - current$summary$theta))
    cv_diff <- max(abs(base$summary$cv_error - current$summary$cv_error))
    lambda_diff <- abs(base$summary$lambda - current$summary$lambda)
    selected_same <- identical(base$summary$selected, current$summary$selected)
    list(max_abs_diff = max(theta_diff, lambda_diff), cv_error_max_abs_diff = cv_diff, selected_identical = selected_same)
  } else {
    mean_diff <- max(abs(base$summary$mean - current$summary$mean))
    selected_same <- identical(base$summary$selected, current$summary$selected)
    list(max_abs_diff = mean_diff, cv_error_max_abs_diff = NA_real_, selected_identical = selected_same)
  }
}

run_impl <- function(label, src) {
  lib <- file.path(tempdir(), paste0("nonprobsvy-lib-", label, "-", timestamp))
  install_pkg(src, lib, label)
  out <- tempfile(fileext = ".rds")
  run_cmd(file.path(R.home("bin"), "Rscript"), c(bench_worker, lib, out, as.character(reps), as.character(quick)))
  readRDS(out)
}

message("Preparing and installing baseline: ", base_ref)
base_src <- prepare_source("baseline", base_ref)
base_rows <- run_impl("baseline", base_src)

message("Preparing and installing current worktree")
current_src <- prepare_source("current")
current_rows <- run_impl("current", current_src)

base_summary <- summarize_one("baseline", base_rows)
current_summary <- summarize_one("current", current_rows)
merged <- merge(
  base_summary,
  current_summary,
  by = c("scenario_id", "estimator", "link", "size", "pop_totals", "est_method"),
  suffixes = c("_baseline", "_current"),
  sort = FALSE
)
merged$speedup <- merged$median_sec_baseline / merged$median_sec_current

num <- lapply(seq_along(base_rows), function(i) numeric_summary(base_rows[[i]], current_rows[[i]]))
merged$max_abs_diff <- vapply(num, `[[`, numeric(1), "max_abs_diff")
merged$cv_error_max_abs_diff <- vapply(num, `[[`, numeric(1), "cv_error_max_abs_diff")
merged$selected_identical <- vapply(num, `[[`, logical(1), "selected_identical")

direct <- merged[merged$estimator == "direct_cv", , drop = FALSE]
e2e <- merged[merged$estimator == "end_to_end", , drop = FALSE]
median_direct_speedup <- median(direct$speedup)
numeric_pass <- all(merged$max_abs_diff <= 1e-6) && all(merged$selected_identical)
end_to_end_regression_free <- all(e2e$speedup >= 0.95)
gate_pass <- median_direct_speedup >= threshold && numeric_pass && end_to_end_regression_free

write.csv(merged, results_csv, row.names = FALSE)
file.copy(results_csv, latest_csv, overwrite = TRUE)

report <- c(
  "# Rcpp CV Benchmark Report",
  "",
  sprintf("- Timestamp: `%s`", timestamp),
  sprintf("- Baseline ref: `%s`", base_ref),
  sprintf("- Current ref: worktree at `%s`", repo_root),
  sprintf("- Repetitions per scenario: `%s`", reps),
  sprintf("- Quick mode: `%s`", quick),
  sprintf("- Required median direct-CV speedup: `%.2fx`", threshold),
  "",
  "## Gate",
  "",
  sprintf("- Median direct-CV speedup: `%.3fx`", median_direct_speedup),
  sprintf("- Numerical equivalence pass (`max_abs_diff <= 1e-6` for fitted outputs and identical selected sets): `%s`", numeric_pass),
  sprintf("- End-to-end regression-free pass (`speedup >= 0.95` for representative workloads): `%s`", end_to_end_regression_free),
  sprintf("- Overall gate: `%s`", if (gate_pass) "PASS" else "FAIL"),
  "",
  "## Scenario Results",
  "",
  paste(
    sprintf(
      "- `%s`: baseline `%.4fs`, current `%.4fs`, speedup `%.3fx`, max fitted diff `%g`, CV-error diff `%s`, selected identical `%s`",
      merged$scenario_id,
      merged$median_sec_baseline,
      merged$median_sec_current,
      merged$speedup,
      merged$max_abs_diff,
      ifelse(is.na(merged$cv_error_max_abs_diff), "NA", format(merged$cv_error_max_abs_diff, scientific = TRUE)),
      merged$selected_identical
    ),
    collapse = "\n"
  ),
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
message("Overall gate: ", if (gate_pass) "PASS" else "FAIL")

if (!gate_pass) quit(status = 2L)
