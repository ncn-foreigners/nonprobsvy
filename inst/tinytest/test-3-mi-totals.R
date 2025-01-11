# These tests are only supposed to be run on developer's machine and
# package GitHub page not on CRAN (they take too long)

if (isTRUE(tolower(Sys.getenv("TEST_NONPROBSVY_MULTICORE_DEVELOPER")) == "true")) {

# Load necessary libraries
library(sampling)
library(survey)

# Generation of data ----------------------------------------------------------------------
set.seed(2024)
source("test-1-generate-data.R")  # Adjust the path as needed

# Check Mass Imputation (MI) Estimator -----------------------------------------------------

## Linear case ------------------------------------------------------------------------------

#### Correctly specified variables ---------------------------------------------------------

##### One target variable ------------------------------------------------------------------

# For Y_11
expect_silent(
  y11_mi_corr_one <- nonprob(
    outcome = Y_11 ~ X1 + X2 + X3 + X4,
    data = sample_B1,
    pop_totals = X_totals[1:5],
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
)

expect_equal(y11_mi_corr_one$output$mean, 2.054878, tolerance = 0.0001)  # Adjusted expected mean
expect_equal(y11_mi_corr_one$output$SE, 0.058286, tolerance = 0.0001)    # Adjusted expected SE
expect_true(y11_mi_corr_one$confidence_interval$lower_bound < mean(Y_11) &
              y11_mi_corr_one$confidence_interval$upper_bound > mean(Y_11))  # Confidence interval test

# For Y_12
expect_silent(
  y12_mi_corr_one <- nonprob(
    outcome = Y_12 ~ X1 + X2 + X3 + X4,
    data = sample_B1,
    pop_totals = X_totals[1:5],
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
)

expect_equal(y12_mi_corr_one$output$mean, 7.465879, tolerance = 0.0001)
expect_equal(y12_mi_corr_one$output$SE, 0.2523682, tolerance = 0.0001)
# expect_true(y12_mi_corr_one$confidence_interval$lower_bound < mean(Y_12) &
#               y12_mi_corr_one$confidence_interval$upper_bound > mean(Y_12))

# For Y_21 (Binary outcome)
expect_silent(
  y21_mi_corr_one <- nonprob(
    outcome = Y_21 ~ X1 + X2 + X3 + X4,
    data = sample_B1,
    pop_totals = X_totals[1:5],
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)

expect_equal(y21_mi_corr_one$output$mean, 0.6706463, tolerance = 0.0001)
expect_equal(y21_mi_corr_one$output$SE, 0.01796104, tolerance = 0.0001)
expect_true(y21_mi_corr_one$confidence_interval$lower_bound < mean(Y_21) &
              y21_mi_corr_one$confidence_interval$upper_bound > mean(Y_21))

# For Y_22 (Binary outcome)
expect_silent(
  y22_mi_corr_one <- nonprob(
    outcome = Y_22 ~ X1 + X2 + X3 + X4,
    data = sample_B1,
    pop_totals = X_totals[1:5],
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)

expect_equal(y22_mi_corr_one$output$mean, 0.653369, tolerance = 0.0001)
expect_equal(y22_mi_corr_one$output$SE, 0.01654733, tolerance = 0.0001)
expect_true(y22_mi_corr_one$confidence_interval$lower_bound < mean(Y_22) &
              y22_mi_corr_one$confidence_interval$upper_bound > mean(Y_22))

##### All target variables -----------------------------------------------------------------

#### All X variables -----------------------------------------------------------------------

##### One target variable ------------------------------------------------------------------

# For Y_11 with all X variables
expect_silent(
  y11_mi_corr_all <- nonprob(
    outcome = as.formula(paste('Y_11', as.character(X_formula))),
    data = sample_B1,  # Include all X variables
    pop_totals = X_totals,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
)

expect_equal(y11_mi_corr_all$output$mean, 2.002616, tolerance = 0.0001)
expect_equal(y11_mi_corr_all$output$SE, 0.03326564, tolerance = 0.0001)
expect_true(y11_mi_corr_all$confidence_interval$lower_bound < mean(Y_11) &
              y11_mi_corr_all$confidence_interval$upper_bound > mean(Y_11))

# For Y_12 with all X variables
expect_silent(
  y12_mi_corr_all <- nonprob(
    outcome = as.formula(paste('Y_12', as.character(X_formula))),
    data = sample_B1,
    pop_totals = X_totals,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
)

expect_equal(y12_mi_corr_all$output$mean, 7.425809, tolerance = 0.0001)
expect_equal(y12_mi_corr_all$output$SE, 0.2454465, tolerance = 0.0001)
# TODO
# expect_true(y12_mi_corr_all$confidence_interval$lower_bound < mean(Y_12) &
#               y12_mi_corr_all$confidence_interval$upper_bound > mean(Y_12))

# For Y_21 with all X variables
expect_silent(
  y21_mi_corr_all <- nonprob(
    outcome = as.formula(paste('Y_21', as.character(X_formula))),
    data = sample_B1,
    pop_totals = X_totals,
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)

expect_equal(y21_mi_corr_all$output$mean, 0.7112867, tolerance = 0.0001)
expect_equal(y21_mi_corr_all$output$SE, 0.01992489, tolerance = 0.0001)
# TODO
# expect_true(y21_mi_corr_all$confidence_interval$lower_bound < mean(Y_21) &
#               y21_mi_corr_all$confidence_interval$upper_bound > mean(Y_21))

# For Y_22 with all X variables
expect_silent(
  y22_mi_corr_all <- nonprob(
    outcome = as.formula(paste('Y_22', as.character(X_formula))),
    data = sample_B1,
    pop_totals = X_totals,
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)

expect_equal(y22_mi_corr_all$output$mean, 0.7836929, tolerance = 0.0001)
expect_equal(y22_mi_corr_all$output$SE, 0.01921207, tolerance = 0.0001)
# TODO
# expect_true(y22_mi_corr_all$confidence_interval$lower_bound < mean(Y_22) &
#               y22_mi_corr_all$confidence_interval$upper_bound > mean(Y_22))

##### All target variables -----------------------------------------------------------------

#### Variable Selection --------------------------------------------------------------------

##### One target variable ------------------------------------------------------------------

# For Y_11 with SCAD penalty
expect_silent(
  y11_mi_scad <- nonprob(
    outcome = Y_11 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    data = sample_B1,
    pop_totals = X_totals[1:11],
    method_outcome = "glm",
    family_outcome = "gaussian",
    control_outcome = controlOut(penalty = "SCAD", nfolds = 5),
    control_inference = controlInf(vars_selection = TRUE)
  )
)

expect_equal(y11_mi_scad$output$mean, 2.00172, tolerance = 0.0001)
expect_equal(y11_mi_scad$output$SE, 0.03371772, tolerance = 0.0001)
expect_true(y11_mi_scad$confidence_interval$lower_bound < mean(Y_11) &
              y11_mi_scad$confidence_interval$upper_bound > mean(Y_11))
expect_true(NROW(y11_mi_scad$outcome$coefficients) < 11)  # Number of selected variables

# For Y_12 with SCAD penalty
expect_silent(
  y12_mi_scad <- nonprob(
    outcome = Y_12 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    data = sample_B1,
    pop_totals = X_totals[1:11],
    method_outcome = "glm",
    family_outcome = "gaussian",
    control_outcome = controlOut(penalty = "SCAD", nfolds = 5),
    control_inference = controlInf(vars_selection = TRUE)
  )
)

expect_equal(y12_mi_scad$output$mean, 7.48752, tolerance = 0.0001)
expect_equal(y12_mi_scad$output$SE, 0.2202497, tolerance = 0.0001)
# TODO
# expect_true(y12_mi_scad$confidence_interval$lower_bound < mean(Y_12) &
#               y12_mi_scad$confidence_interval$upper_bound > mean(Y_12))
expect_true(NROW(y12_mi_scad$outcome$coefficients) < 11)

# For Y_21 with SCAD penalty
expect_silent(
  y21_mi_scad <- nonprob(
    outcome = Y_21 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    data = sample_B1,
    pop_totals = X_totals[1:11],
    method_outcome = "glm",
    family_outcome = "binomial",
    control_outcome = controlOut(penalty = "SCAD", nfolds = 5),
    control_inference = controlInf(vars_selection = TRUE)
  )
)

expect_equal(y21_mi_scad$output$mean, 0.7158588, tolerance = 0.0001)
expect_equal(y21_mi_scad$output$SE, 0.0141177, tolerance = 0.0001)
# TODO
# expect_true(y21_mi_scad$confidence_interval$lower_bound < mean(Y_21) &
#               y21_mi_scad$confidence_interval$upper_bound > mean(Y_21))
expect_true(NROW(y21_mi_scad$outcome$coefficients) < 11)

# For Y_22 with SCAD penalty
expect_silent(
  y22_mi_scad <- nonprob(
    outcome = Y_22 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    data = sample_B1,
    pop_totals = X_totals[1:11],
    method_outcome = "glm",
    family_outcome = "binomial",
    control_outcome = controlOut(penalty = "SCAD", nfolds = 5),
    control_inference = controlInf(vars_selection = TRUE)
  )
)

expect_equal(y22_mi_scad$output$mean, 0.7870636, tolerance = 0.0001)
expect_equal(y22_mi_scad$output$SE, 0.01550724, tolerance = 0.0001)
# TODO
# expect_true(y22_mi_scad$confidence_interval$lower_bound < mean(Y_22) &
#               y22_mi_scad$confidence_interval$upper_bound > mean(Y_22))
expect_true(NROW(y22_mi_scad$outcome$coefficients) < 11)
}
