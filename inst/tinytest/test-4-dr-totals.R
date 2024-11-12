# Load necessary libraries
library(sampling)
library(survey)

# Generation of data ----------------------------------------------------------------------
set.seed(2024)
source("test-1-generate-data.R")  # Adjust the path as needed

# Check Doubly Robust (DR) Estimator -------------------------------------------------------

## Linear case -----------------------------------------------------------------------------

#### Correctly specified variables ---------------------------------------------------------

##### One target variable ------------------------------------------------------------------

# For Y_11
expect_silent(
  y11_dr_corr_one <- nonprob(
    selection = ~ X1 + X2 + X3 + X4,
    outcome = Y_11 ~ X1 + X2 + X3 + X4,
    data = sample_B1,
    pop_totals = X_totals[1:5],
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
)

expect_equal(y11_dr_corr_one$output$mean, 2.17757, tolerance = 0.0001)  # Adjusted expected mean
expect_equal(y11_dr_corr_one$output$SE, 0.09009875, tolerance = 0.0001)    # Adjusted expected SE
# TODO
# expect_true(y11_dr_corr_one$confidence_interval$lower_bound < mean(Y_11) &
#               y11_dr_corr_one$confidence_interval$upper_bound > mean(Y_11))  # Confidence interval test

# For Y_12
expect_silent(
  y12_dr_corr_one <- nonprob(
    selection = ~ X1 + X2 + X3 + X4,
    outcome = Y_12 ~ X1 + X2 + X3 + X4,
    data = sample_B1,
    pop_totals = X_totals[1:5],
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
)

expect_equal(y12_dr_corr_one$output$mean, 7.168049, tolerance = 0.0001)
expect_equal(y12_dr_corr_one$output$SE, 0.5737179, tolerance = 0.0001)
expect_true(y12_dr_corr_one$confidence_interval$lower_bound < mean(Y_12) &
              y12_dr_corr_one$confidence_interval$upper_bound > mean(Y_12))

# For Y_21 (Binary outcome)
expect_silent(
  y21_dr_corr_one <- nonprob(
    selection = ~ X1 + X2 + X3 + X4,
    outcome = Y_21 ~ X1 + X2 + X3 + X4,
    data = sample_B1,
    pop_totals = X_totals[1:5],
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)

expect_equal(y21_dr_corr_one$output$mean, 0.6969967, tolerance = 0.0001)
expect_equal(y21_dr_corr_one$output$SE, 0.02966577, tolerance = 0.0001)
expect_true(y21_dr_corr_one$confidence_interval$lower_bound < mean(Y_21) &
              y21_dr_corr_one$confidence_interval$upper_bound > mean(Y_21))

# For Y_22 (Binary outcome)
expect_silent(
  y22_dr_corr_one <- nonprob(
    selection = ~ X1 + X2 + X3 + X4,
    outcome = Y_22 ~ X1 + X2 + X3 + X4,
    data = sample_B1,
    pop_totals = X_totals[1:5],
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)

expect_equal(y22_dr_corr_one$output$mean, 0.62005, tolerance = 0.0001)
expect_equal(y22_dr_corr_one$output$SE, 0.05059186, tolerance = 0.0001)
expect_true(y22_dr_corr_one$confidence_interval$lower_bound < mean(Y_22) &
              y22_dr_corr_one$confidence_interval$upper_bound > mean(Y_22))

##### All target variables -----------------------------------------------------------------

#### All X variables -----------------------------------------------------------------------

##### One target variable ------------------------------------------------------------------

# For Y_11 with all X variables
expect_silent(
  y11_dr_corr_all <- nonprob(
    selection = X_formula,
    outcome = as.formula(paste('Y_11 ~ ', as.character(X_formula)[2])),
    data = sample_B1,
    pop_totals = X_totals,
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
)

expect_equal(y11_dr_corr_all$output$mean, 2.005841, tolerance = 0.0001)
expect_equal(y11_dr_corr_all$output$SE, 0.04629972, tolerance = 0.0001)
expect_true(y11_dr_corr_all$confidence_interval$lower_bound < mean(Y_11) &
              y11_dr_corr_all$confidence_interval$upper_bound > mean(Y_11))

# For Y_12 with all X variables
expect_silent(
  y12_dr_corr_all <- nonprob(
    selection = X_formula,
    outcome = as.formula(paste('Y_12 ~ ', as.character(X_formula)[2])),
    data = sample_B1,
    pop_totals = X_totals,
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
)

expect_equal(y12_dr_corr_all$output$mean, 6.681308, tolerance = 0.0001)
expect_equal(y12_dr_corr_all$output$SE, 0.3886612, tolerance = 0.0001)
expect_true(y12_dr_corr_all$confidence_interval$lower_bound < mean(Y_12) &
              y12_dr_corr_all$confidence_interval$upper_bound > mean(Y_12))

# For Y_21 with all X variables
expect_silent(
  y21_dr_corr_all <- nonprob(
    selection = X_formula,
    outcome = as.formula(paste('Y_21 ~ ', as.character(X_formula)[2])),
    data = sample_B1,
    pop_totals = X_totals,
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)

expect_equal(y21_dr_corr_all$output$mean, 0.7233015, tolerance = 0.0001)
expect_equal(y21_dr_corr_all$output$SE, 0.0192478, tolerance = 0.0001)
# TODO
# expect_true(y21_dr_corr_all$confidence_interval$lower_bound < mean(Y_21) &
#               y21_dr_corr_all$confidence_interval$upper_bound > mean(Y_21))

# For Y_22 with all X variables
expect_silent(
  y22_dr_corr_all <- nonprob(
    selection = X_formula,
    outcome = as.formula(paste('Y_22 ~ ', as.character(X_formula)[2])),
    data = sample_B1,
    pop_totals = X_totals,
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "binomial"
  )
)

expect_equal(y22_dr_corr_all$output$mean, 0.7668045, tolerance = 0.0001)
expect_equal(y22_dr_corr_all$output$SE, 0.02748494, tolerance = 0.0001)
# TODO
# expect_true(y22_dr_corr_all$confidence_interval$lower_bound < mean(Y_22) &
#               y22_dr_corr_all$confidence_interval$upper_bound > mean(Y_22))

##### All target variables -----------------------------------------------------------------


#### Variable Selection --------------------------------------------------------------------

##### One target variable ------------------------------------------------------------------

# For Y_11 with SCAD penalty
expect_silent(
  y11_dr_scad <- nonprob(
    selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    outcome = Y_11 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    data = sample_B1,
    pop_totals = X_totals[1:11],
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "gaussian",
    control_selection = controlSel(penalty = "SCAD", nfolds = 5),
    control_outcome = controlOut(penalty = "SCAD", nfolds = 5),
    control_inference = controlInf(vars_selection = TRUE)
  )
)

expect_equal(y11_dr_scad$output$mean, 1.975596, tolerance = 0.0001)
expect_equal(y11_dr_scad$output$SE, 0.03277463, tolerance = 0.0001)
expect_true(y11_dr_scad$confidence_interval$lower_bound < mean(Y_11) &
              y11_dr_scad$confidence_interval$upper_bound > mean(Y_11))
expect_true(NROW(y11_dr_scad$selection$coefficients) < 11)
expect_true(NROW(y11_dr_scad$outcome$coefficients) < 11)

# For Y_12 with SCAD penalty
expect_silent(
  y12_dr_scad <- nonprob(
    selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    outcome = Y_12 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    data = sample_B1,
    pop_totals = X_totals[1:11],
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "gaussian",
    control_selection = controlSel(penalty = "SCAD", nfolds = 5),
    control_outcome = controlOut(penalty = "SCAD", nfolds = 5),
    control_inference = controlInf(vars_selection = TRUE)
  )
)

expect_equal(y12_dr_scad$output$mean, 6.730679, tolerance = 0.0001)
expect_equal(y12_dr_scad$output$SE, 0.5256198, tolerance = 0.0001)
expect_true(y12_dr_scad$confidence_interval$lower_bound < mean(Y_12) &
              y12_dr_scad$confidence_interval$upper_bound > mean(Y_12))
expect_true(NROW(y12_dr_scad$selection$coefficients) < 11)
expect_true(NROW(y12_dr_scad$outcome$coefficients) < 11)

# For Y_21 with SCAD penalty
expect_silent(
  y21_dr_scad <- nonprob(
    selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    outcome = Y_21 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    data = sample_B1,
    pop_totals = X_totals[1:11],
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "binomial",
    control_selection = controlSel(penalty = "SCAD", nfolds = 5),
    control_outcome = controlOut(penalty = "SCAD", nfolds = 5),
    control_inference = controlInf(vars_selection = TRUE)
  )
)

expect_equal(y21_dr_scad$output$mean, 0.7286048, tolerance = 0.0001)
expect_equal(y21_dr_scad$output$SE, 0.01451463, tolerance = 0.0001)
# TODO
# expect_true(y21_dr_scad$confidence_interval$lower_bound < mean(Y_21) &
#               y21_dr_scad$confidence_interval$upper_bound > mean(Y_21))
expect_true(NROW(y21_dr_scad$selection$coefficients) < 11)
expect_true(NROW(y21_dr_scad$outcome$coefficients) < 11)

# For Y_22 with SCAD penalty
expect_silent(
  y22_dr_scad <- nonprob(
    selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    outcome = Y_22 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
    data = sample_B1,
    pop_totals = X_totals[1:11],
    method_selection = "logit",
    method_outcome = "glm",
    family_outcome = "binomial",
    control_selection = controlSel(penalty = "SCAD", nfolds = 5),
    control_outcome = controlOut(penalty = "SCAD", nfolds = 5),
    control_inference = controlInf(vars_selection = TRUE)
  )
)

expect_equal(y22_dr_scad$output$mean, 0.7644814, tolerance = 0.0001)
expect_equal(y22_dr_scad$output$SE, 0.01492712, tolerance = 0.0001)
# TODO
# expect_true(y22_dr_scad$confidence_interval$lower_bound < mean(Y_22) &
#               y22_dr_scad$confidence_interval$upper_bound > mean(Y_22))
expect_true(NROW(y22_dr_scad$selection$coefficients) < 11)
expect_true(NROW(y22_dr_scad$outcome$coefficients) < 11)
