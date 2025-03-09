# no print doccumentation


# general printing --------------------------------------------------------

#' @method print nonprob
#' @exportS3Method
print.nonprob <- function(x, digits=4,...) {

  cat(sprintf("A %s object\n", class(x)[1L]))
  cat(sprintf(" - estimator type: %s\n", switch(x$estimator,
                                                "mi"= "mass imputation",
                                                "ipw" = "inverse probability weighting",
                                                "dr" = "doubly robust")))

  ## here we should have more information when the doubly robust estimator is applied

  cat(sprintf(" - method: %s\n", x$estimator_method))
  cat(sprintf(" - auxiliary variables source: %s\n", ifelse(!is.null(x$svydesign), "survey", "population")))
  cat(sprintf(" - vars selection: %s\n", tolower(x$control$control_inference$vars_selection)))

  if (x$control$control_inference$var_method == "analytic") {
    cat(sprintf(" - variance estimator: analytic\n"))
  } else {
    cat(sprintf(" - variance estimator: bootstrap (with %s samples, multicore: %s)\n",
                x$control$control_inference$num_boot,
                tolower(x$control$control_inference$cores > 1)))
  }

  cat(sprintf(" - population size fixed: %s\n", tolower(x$pop_size_fixed)))

  if (NROW(x$output) == 1) {
    cat(sprintf(" - naive (uncorrected) estimator: %.*f\n", digits, as.numeric(mean(x$y[[1]]))))
    cat(sprintf(
      " - selected estimator: %.*f (se=%.*f, ci=(%.*f, %.*f))\n",
      digits, as.numeric(x$output["mean"]),
      digits, as.numeric(x$output["SE"]),
      digits, as.numeric(x$confidence_interval["lower_bound"]),
      digits, as.numeric(x$confidence_interval["upper_bound"])))
  } else if (NROW(x$output) <= 5) {
    cat(sprintf(" - naive (uncorrected) estimators:\n"))
    for (v in 1:NROW(x$output)) {
          cat(sprintf("   - variable %s: %.*f\n",
              rownames(x$output)[v],
              digits,
              as.numeric(mean(x$y[[v]]))
              ))
    }

    cat(" - selected estimators:\n")
    for (v in 1:NROW(x$output)) {
      cat(sprintf(
        "   - variable %s: %.*f (se=%.*f, ci=(%.*f, %.*f))\n",
        rownames(x$output)[v],
        digits, as.numeric(x$output[v, "mean"]),
        digits, as.numeric(x$output[v, "SE"]),
        digits, as.numeric(x$confidence_interval[v, "lower_bound"]),
        digits, as.numeric(x$confidence_interval[v,"upper_bound"])))
    }
  }
  else {
    cat(" - too many variables (more than 5) to summary. See details in the `output` and `confidence_interval` elements.\n")
  }

  #cat("Estimated population mean with overall std.err and confidence interval:\n\n")
  #print(cbind(mean = x$output$mean, SE = x$output$SE, x$confidence_interval))
  invisible(x)
}
#' @title Print method for the `nonprob_summary` object
#'
#' @description
#' Print method for the `nonprob_summary` object which allows for specification
#' what should be printed or not.
#'
#' @param x a `nonprob` object
#' @param resid whether distribution of residuals should be printed (default is `TRUE`)
#' @param pred whether distribution of predictions should be printed (default is `TRUE`)
#' @param digits number of digits to be printed (default 4)
#' @param ... further parameters passed to the print method (currently not supported)
#'
#' @method print nonprob_summary
#' @exportS3Method
print.nonprob_summary <- function(x,
                                  resid = TRUE,
                                  pred = TRUE,
                                  digits = 4,
                                  ...) {
  cat(sprintf("A %s object\n", class(x)[1L]))
  if (!is.null(x$call)) {
    cat(" - call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n",
        sep = ""
    )
  }
  cat(sprintf(" - estimator type: %s\n", switch(x$estimator,
                                                "mi"= "mass imputation",
                                                "ipw" = "inverse probability weighting",
                                                "dr" = "doubly robust")))

  cat(" - nonprob sample size: ", x$nonprob_size, " (", round(x$nonprob_size/x$pop_size*100,1),"%)\n", sep = "")

  if (x$prob_size > 0) {
    cat(" - prob sample size: ", x$prob_size, " (", round(x$prob_size/x$pop_size*100,1),"%)\n", sep = "")
  } else {
    cat(" - prob sample size: none; population totals were provided\n")
  }

  cat(" - population size: ", as.integer(x$pop_size), " (fixed: ", tolower(x$pop_size_fixed), ")\n", sep = "")

  model_info <- switch(x$estimator,
                       "mi" = '"outcome"',
                       "ipw" = '"selection"',
                       "dr" = '"outcome" and "selection"')
  cat(sprintf(
    " - detailed information about models are stored in list element(s): %s\n",
    model_info))

  if (!is.null(x$ipw_weights)) {
    cat("----------------------------------------------------------------\n")
      cat(" - sum of IPW weights:", sum(x$ipw_weights_total), "\n")
      cat(" - distribution of IPW weights (nonprob sample):\n")
      cat(sprintf(
        "   - min: %.*f; mean: %.*f; median: %.*f; max: %.*f\n",
        digits, as.numeric(x$ipw_weights["Min."]),
        digits, as.numeric(x$ipw_weights["Mean"]),
        digits, as.numeric(x$ipw_weights["Median"]),
        digits, as.numeric(x$ipw_weights["Max."])))

      cat(" - distribution of IPW probabilities (nonprob sample):\n")
      cat(sprintf(
        "   - min: %.*f; mean: %.*f; median: %.*f; max: %.*f\n",
        digits, as.numeric(x$ps_scores_nonprob["Min."]),
        digits, as.numeric(x$ps_scores_nonprob["Mean"]),
        digits, as.numeric(x$ps_scores_nonprob["Median"]),
        digits, as.numeric(x$ps_scores_nonprob["Max."])))

      if (!x$no_prob) {
        cat(" - distribution of IPW probabilities (prob sample):\n")
        cat(sprintf(
          "   - min: %.*f; mean: %.*f; median: %.*f; max: %.*f\n",
          digits, as.numeric(x$ps_scores_prob["Min."]),
          digits, as.numeric(x$ps_scores_prob["Mean"]),
          digits, as.numeric(x$ps_scores_prob["Median"]),
          digits, as.numeric(x$ps_scores_prob["Max."])))
      }
      if (is.null(x$ys_resid)) {
        cat("----------------------------------------------------------------\n")
      }
  }

  if (!is.null(x$ys_resid)) {
    cat("----------------------------------------------------------------\n")
    if (resid) {
      if (length(x$ys_resid) <= 5) {
        cat(" - distribution of outcome residuals:\n")
        for (y in 1:length(x$ys_resid)) {
          cat(sprintf(
            "   - %s: min: %.*f; mean: %.*f; median: %.*f; max: %.*f\n",
            names(x$ys_resid)[y],
            digits, as.numeric(x$ys_resid[[y]]["Min."]),
            digits, as.numeric(x$ys_resid[[y]]["Mean"]),
            digits, as.numeric(x$ys_resid[[y]]["Median"]),
            digits, as.numeric(x$ys_resid[[y]]["Max."])))
        }
      } else {
        cat(" - distribution of outcome residuals:\n")
        cat("   too many variables (more than 5). Details are stored in the `ys_resid` element.")
      }

    }
    if (pred) {
      if (length(x$ys_nons_pred) <= 5) {
        cat(" - distribution of outcome predictions (nonprob sample):\n")
        for (y in 1:length(x$ys_nons_pred)) {
          cat(sprintf(
            "   - %s: min: %.*f; mean: %.*f; median: %.*f; max: %.*f\n",
            names(x$ys_nons_pred)[y],
            digits, as.numeric(x$ys_nons_pred[[y]]["Min."]),
            digits, as.numeric(x$ys_nons_pred[[y]]["Mean"]),
            digits, as.numeric(x$ys_nons_pred[[y]]["Median"]),
            digits, as.numeric(x$ys_nons_pred[[y]]["Max."])))
        }

        if (!x$no_prob) {
          cat(" - distribution of outcome predictions (prob sample):\n")
          for (y in 1:length(x$ys_rand_pred)) {
            cat(sprintf(
              "   - %s: min: %.*f; mean: %.*f; median: %.*f; max: %.*f\n",
              names(x$ys_rand_pred)[y],
              digits, as.numeric(x$ys_rand_pred[[y]]["Min."]),
              digits, as.numeric(x$ys_rand_pred[[y]]["Mean"]),
              digits, as.numeric(x$ys_rand_pred[[y]]["Median"]),
              digits, as.numeric(x$ys_rand_pred[[y]]["Max."])))
          }
        }

      } else {
        cat(" - distribution of outcome predictions:\n")
        cat("   too many variables (more than 5). Details are stored in the `ys_nons_pred` and `ys_rand_pred` elements.")
      }
    }
    cat("----------------------------------------------------------------\n")
    }

  invisible(x)
}



# print for specific classes ------------------------------------------------

#' @method print nonprob_method
#' @exportS3Method
print.nonprob_method <- function(x, ...) {

  if (x$model == "ps") {
    cat(sprintf("Propensity score model with %s link", x$link))
  } else {
    cat(sprintf("Mass imputation model (%s approach). Estimated mean: %s (se: %s)",
                  toupper(x$model),
                  sprintf("%.4f", x$y_mi_hat),
                  sprintf("%.4f", sqrt(x$var_prob+x$var_nonprob))))
  }

  invisible(x)
}


