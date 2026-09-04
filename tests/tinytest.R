if (requireNamespace("tinytest", quietly = TRUE)) {
  full_tests <- identical(tolower(Sys.getenv("NONPROBSVY_FULL_TESTS")), "true")

  if (full_tests) {
    tinytest::test_package("nonprobsvy")
  } else {
    message(
      "Running the CRAN-fast tinytest subset. ",
      "Set NONPROBSVY_FULL_TESTS=true to run the full suite."
    )
    fast_tests <- c(
      "test_check_balance",
      "test_controls",
      "test_ipw_stability",
      "test_ipw_mle_variance_formulas",
      "test_make_outcomes",
      "test_methods",
      "test_mi",
      "test_mi_known_n",
      "test_mi_methods",
      "test_nn_cache",
      "test_nonprob",
      "test_nonprob_unit_level",
      "test_pmm_k_choice",
      "test_rcpp_cv"
    )
    tinytest::test_package(
      "nonprobsvy",
      pattern = paste0("^(", paste(fast_tests, collapse = "|"), ")\\.R$")
    )
  }
}
