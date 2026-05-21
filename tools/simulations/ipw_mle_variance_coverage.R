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

reps <- as.integer(parse_arg("reps", "500"))
seed <- as.integer(parse_arg("seed", "20260509"))
links <- strsplit(parse_arg("links", "logit,probit,cloglog"), ",", fixed = TRUE)[[1]]
links <- trimws(links)
known_n <- identical(parse_arg("known-n", "true"), "true")

if (anyNA(reps) || reps < 1) {
  stop("`--reps` must be a positive integer.", call. = FALSE)
}

if (any(!links %in% c("logit", "probit", "cloglog"))) {
  stop("`--links` must contain only logit, probit, cloglog.", call. = FALSE)
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
    logit = -1.45 + 0.55 * x1 - 0.35 * x2,
    probit = -0.85 + 0.35 * x1 - 0.20 * x2,
    cloglog = -1.70 + 0.45 * x1 - 0.30 * x2
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

run_one <- function(pop, link, known_n) {
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
      pop_size = if (known_n) nrow(pop) else NULL,
      method_selection = link,
      control_selection = control_sel(
        est_method = "mle",
        optimizer = "optim",
        optim_method = "BFGS",
        maxit = 1000
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

summarise_results <- function(results) {
  do.call(rbind, lapply(split(results, results$link), function(x) {
    err <- x$estimate - x$true_mean
    data.frame(
      link = x$link[[1]],
      reps = nrow(x),
      failures = NA_integer_,
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
}

main <- function() {
  pop <- make_population(seed)
  set.seed(seed + 1L)

  all_results <- list()
  attempted <- setNames(integer(length(links)), links)
  for (link in links) {
    for (r in seq_len(reps)) {
      attempted[[link]] <- attempted[[link]] + 1L
      result <- run_one(pop, link, known_n = known_n)
      if (!is.null(result)) {
        all_results[[length(all_results) + 1L]] <- result
      }
      if (r %% 50L == 0L) {
        message(sprintf("completed %s: %d/%d", link, r, reps))
      }
    }
  }

  results <- do.call(rbind, all_results)
  summary <- summarise_results(results)
  summary$failures <- attempted[summary$link] - summary$reps

  cat("Simulation design: fixed finite population; nonprobability Bernoulli sample; ")
  cat("stratified one-stage cluster probability sample with strata, PSU weights, and fpc.\n")
  cat(sprintf("Population N: %d; true mean: %.6f; known_n: %s\n\n",
              nrow(pop), mean(pop$y), known_n))
  print(summary, row.names = FALSE, digits = 4)

  coverage_ok <- summary$coverage >= 0.93 & summary$coverage <= 0.97
  se_ok <- summary$se_ratio >= 0.90 & summary$se_ratio <= 1.10

  if (all(coverage_ok & se_ok)) {
    quit(status = 0)
  }

  quit(status = 1)
}

if (sys.nframe() == 0L) {
  main()
}
