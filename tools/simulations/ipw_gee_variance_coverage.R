#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(survey)
})

parse_arg <- function(name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[1]])
}

reps <- as.integer(parse_arg("reps", "300"))
seed <- as.integer(parse_arg("seed", "20260510"))
links <- strsplit(parse_arg("links", "logit,probit,cloglog"), ",", fixed = TRUE)[[1]]
links <- trimws(links)
h_funs <- as.integer(strsplit(parse_arg("h", "1,2"), ",", fixed = TRUE)[[1]])

if (anyNA(reps) || reps < 1) {
  stop("`--reps` must be a positive integer.", call. = FALSE)
}

if (any(!links %in% c("logit", "probit", "cloglog"))) {
  stop("`--links` must contain only logit, probit, cloglog.", call. = FALSE)
}

if (anyNA(h_funs) || any(!h_funs %in% c(1L, 2L))) {
  stop("`--h` must contain only 1 and/or 2.", call. = FALSE)
}

make_population <- function(seed) {
  set.seed(seed)
  strata_clusters <- c(44L, 52L, 48L, 56L)
  sampled_clusters <- c(11L, 13L, 12L, 14L)
  cluster_size <- 30L
  stratum_effect <- c(-0.45, -0.10, 0.20, 0.55)

  pop <- do.call(rbind, lapply(seq_along(strata_clusters), function(h) {
    psu <- seq_len(strata_clusters[[h]])
    out <- do.call(rbind, lapply(psu, function(g) {
      cluster_effect <- rnorm(1, sd = 0.35)
      x1 <- rnorm(cluster_size, mean = (h - 2.5) / 4 + cluster_effect / 3)
      x2 <- rbinom(cluster_size, 1, plogis(-0.1 + 0.35 * x1 + 0.25 * h))
      y <- 1.4 + 0.7 * x1 - 0.45 * x2 + stratum_effect[[h]] +
        cluster_effect + rnorm(cluster_size, sd = 0.8)
      data.frame(
        stratum = h,
        psu = sprintf("%02d-%03d", h, g),
        x1 = x1,
        x2 = x2,
        y = y
      )
    }))
    out$fpc_psu <- strata_clusters[[h]]
    out$sampled_psu <- sampled_clusters[[h]]
    out$prob_weight <- strata_clusters[[h]] / sampled_clusters[[h]]
    out
  }))

  rownames(pop) <- NULL
  pop$stratum <- factor(pop$stratum)
  pop$unit_id <- seq_len(nrow(pop))
  pop
}

link_eta <- function(link, x1, x2) {
  switch(
    link,
    logit = -1.35 + 0.50 * x1 - 0.30 * x2,
    probit = -0.80 + 0.32 * x1 - 0.18 * x2,
    cloglog = -1.65 + 0.42 * x1 - 0.26 * x2
  )
}

link_inv <- function(link, eta) {
  switch(
    link,
    logit = plogis(eta),
    probit = pnorm(eta),
    cloglog = -expm1(-exp(eta))
  )
}

draw_complex_reference <- function(pop) {
  sampled_psus <- unlist(tapply(pop$psu, pop$stratum, function(psu) {
    psu <- unique(psu)
    h <- pop$stratum[match(psu, pop$psu)]
    m_h <- unique(pop$sampled_psu[pop$stratum == h[[1]]])
    sample(psu, size = m_h, replace = FALSE)
  }), use.names = FALSE)

  prob <- pop[pop$psu %in% sampled_psus, c(
    "unit_id", "stratum", "psu", "x1", "x2", "prob_weight", "fpc_psu"
  )]

  survey::svydesign(
    ids = ~psu,
    strata = ~stratum,
    weights = ~prob_weight,
    fpc = ~fpc_psu,
    data = prob,
    nest = TRUE
  )
}

run_one <- function(pop, link, gee_h_fun) {
  eta <- link_eta(link, pop$x1, pop$x2)
  p <- link_inv(link, eta)
  nons <- pop[stats::rbinom(nrow(pop), 1, p) == 1, c("x1", "x2", "y")]
  if (nrow(nons) < 50) {
    return(NULL)
  }

  svy <- draw_complex_reference(pop)

  fit <- try(
    nonprob(
      data = nons,
      selection = ~x1 + x2,
      target = ~y,
      svydesign = svy,
      pop_size = nrow(pop),
      method_selection = link,
      control_selection = control_sel(
        est_method = "gee",
        gee_h_fun = gee_h_fun,
        maxit = 1000,
        start_type = "mle"
      ),
      control_inference = control_inf(var_method = "analytic"),
      verbose = FALSE
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(NULL)
  }

  est <- fit$output$mean[[1]]
  se <- fit$output$SE[[1]]
  ci <- fit$confidence_interval[1, ]
  true_mean <- mean(pop$y)
  data.frame(
    link = link,
    gee_h_fun = gee_h_fun,
    estimate = est,
    se = se,
    covered = ci$lower_bound <= true_mean && ci$upper_bound >= true_mean,
    true_mean = true_mean,
    nonprob_n = nrow(nons),
    prob_n = fit$prob_size,
    prob_se = fit$SE$prob[[1]],
    nonprob_se = fit$SE$nonprob[[1]]
  )
}

summarise_results <- function(results, attempted) {
  split_key <- interaction(results$link, results$gee_h_fun, drop = TRUE)
  out <- do.call(rbind, lapply(split(results, split_key), function(x) {
    err <- x$estimate - x$true_mean
    key <- paste(x$link[[1]], x$gee_h_fun[[1]], sep = ".")
    data.frame(
      link = x$link[[1]],
      gee_h_fun = x$gee_h_fun[[1]],
      reps = nrow(x),
      failures = unname(attempted[[key]]) - nrow(x),
      bias = mean(err),
      empirical_sd = stats::sd(x$estimate),
      mean_se = mean(x$se),
      se_ratio = mean(x$se) / stats::sd(x$estimate),
      coverage = mean(x$covered),
      mean_nonprob_n = mean(x$nonprob_n),
      mean_prob_n = mean(x$prob_n),
      mean_prob_se = mean(x$prob_se),
      mean_nonprob_se = mean(x$nonprob_se)
    )
  }))
  rownames(out) <- NULL
  out[order(out$link, out$gee_h_fun), ]
}

main <- function() {
  pop <- make_population(seed)
  set.seed(seed + 1L)

  grid <- expand.grid(link = links, gee_h_fun = h_funs, stringsAsFactors = FALSE)
  attempted <- setNames(integer(nrow(grid)), paste(grid$link, grid$gee_h_fun, sep = "."))
  all_results <- list()

  for (i in seq_len(nrow(grid))) {
    link <- grid$link[[i]]
    gee_h_fun <- grid$gee_h_fun[[i]]
    key <- paste(link, gee_h_fun, sep = ".")
    for (r in seq_len(reps)) {
      attempted[[key]] <- attempted[[key]] + 1L
      result <- run_one(pop, link, gee_h_fun)
      if (!is.null(result)) {
        all_results[[length(all_results) + 1L]] <- result
      }
      if (r %% 50L == 0L) {
        message(sprintf("completed %s h=%d: %d/%d", link, gee_h_fun, r, reps))
      }
    }
  }

  results <- do.call(rbind, all_results)
  summary <- summarise_results(results, attempted)

  cat("Simulation design: fixed finite population; nonprobability Bernoulli sample; ")
  cat("stratified one-stage cluster probability sample with strata, PSU weights, and fpc.\n")
  cat(sprintf("Population N: %d; true mean: %.6f\n\n", nrow(pop), mean(pop$y)))
  print(summary, row.names = FALSE, digits = 4)

  coverage_mc_se <- sqrt(0.95 * 0.05 / summary$reps)
  coverage_low <- 0.95 - 2 * coverage_mc_se
  coverage_high <- 0.95 + 2 * coverage_mc_se
  coverage_ok <- summary$coverage >= coverage_low & summary$coverage <= coverage_high

  if (all(coverage_ok)) {
    quit(status = 0)
  }

  quit(status = 1)
}

if (sys.nframe() == 0L) {
  main()
}
