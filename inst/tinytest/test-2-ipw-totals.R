library(sampling)
library(survey)


# generation of data ----------------------------------------------------------------------
set.seed(2024)
source("test-1-generate-data.R") ## path should be changed


# check logit -----------------------------------------------------------------------------
## linear case ----------------------------------------------------------------------------
#### correctly specified variables --------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
expect_silent(
  y11_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                      target = ~ Y_11,
                      data = sample_B1,
                      pop_totals = X_totals[1:5],
                      method_selection = "logit")
)

expect_equal(y11_corr_one$output$mean, 2.17757, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_one$output$SE, 0.1553234, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_one$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_one$confidence_interval$upper_bound > mean(Y_11)) ## conf int


## for y12
expect_silent(
  y12_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                      target = ~ Y_12,
                      data = sample_B1,
                      pop_totals = X_totals[1:5],
                      method_selection = "logit")
)
expect_equal(y12_corr_one$output$mean, 7.168049, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_one$output$SE, 1.01074, tolerance = 0.0001) ## true value for this sim
expect_true(y12_corr_one$confidence_interval$lower_bound < mean(Y_12) &
              y12_corr_one$confidence_interval$upper_bound > mean(Y_12)) ## conf int


## for y21
expect_silent(
  y21_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                      target = ~ Y_21,
                      data = sample_B1,
                      pop_totals = X_totals[1:5],
                      method_selection = "logit")
)

expect_equal(y21_corr_one$output$mean, 0.6675022, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_one$output$SE, 0.05836787, tolerance = 0.0001) ## true value for this sim
expect_true(y21_corr_one$confidence_interval$lower_bound < mean(Y_21) &
              y21_corr_one$confidence_interval$upper_bound > mean(Y_21)) ## conf int

## for y22
expect_silent(
  y22_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                      target = ~ Y_22,
                      data = sample_B1,
                      pop_totals = X_totals[1:5],
                      method_selection = "logit")
)

expect_equal(y22_corr_one$output$mean, 0.6117923, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_one$output$SE, 0.06099473, tolerance = 0.0001) ## true value for this sim
expect_true(y22_corr_one$confidence_interval$lower_bound < mean(Y_22) &
              y22_corr_one$confidence_interval$upper_bound > mean(Y_22)) ## conf int

##### all target variables  ---------------------------------------------------------------

expect_silent(
  y_all_corr <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                        target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                        data = sample_B1,
                        pop_totals = X_totals[1:5],
                        method_selection = "logit")
)

expect_identical(y_all_corr$output$mean,
                 c(y11_corr_one$output$mean, y12_corr_one$output$mean,
                   y21_corr_one$output$mean, y22_corr_one$output$mean))

expect_identical(y_all_corr$output$SE,
                 c(y11_corr_one$output$SE, y12_corr_one$output$SE,
                   y21_corr_one$output$SE, y22_corr_one$output$SE))

expect_identical(y_all_corr$confidence_interval,
                 data.frame(lower_bound = c(y11_corr_one$confidence_interval$lower_bound,
                                            y12_corr_one$confidence_interval$lower_bound,
                                            y21_corr_one$confidence_interval$lower_bound,
                                            y22_corr_one$confidence_interval$lower_bound),
                            upper_bound = c(y11_corr_one$confidence_interval$upper_bound,
                                            y12_corr_one$confidence_interval$upper_bound,
                                            y21_corr_one$confidence_interval$upper_bound,
                                            y22_corr_one$confidence_interval$upper_bound),
                            row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))





#### all X variables variables ------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
expect_silent(
  y11_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_11,
                          data = sample_B1,
                          pop_totals = X_totals,
                          method_selection = "logit")
)

expect_equal(y11_corr_all$output$mean, 2.005841, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_all$output$SE, 0.175094, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_all$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_all$confidence_interval$upper_bound > mean(Y_11)) ## conf int


## for y12
expect_silent(
  y12_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_12,
                          data = sample_B1,
                          pop_totals = X_totals,
                          method_selection = "logit")
)
expect_equal(y12_corr_all$output$mean, 6.681308, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_all$output$SE, 0.7386469, tolerance = 0.0001) ## true value for this sim
expect_true(y12_corr_all$confidence_interval$lower_bound < mean(Y_12) &
              y12_corr_all$confidence_interval$upper_bound > mean(Y_12)) ## conf int


## for y21
expect_silent(
  y21_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_21,
                          data = sample_B1,
                          pop_totals = X_totals,
                          method_selection = "logit")
)

expect_equal(y21_corr_all$output$mean, 0.64798484, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_all$output$SE, 0.061434896, tolerance = 0.0001) ## true value for this sim
expect_true(y21_corr_all$confidence_interval$lower_bound < mean(Y_21) &
              y21_corr_all$confidence_interval$upper_bound > mean(Y_21)) ## conf int

## for y22
expect_silent(
  y22_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_22,
                          data = sample_B1,
                          pop_totals = X_totals,
                          method_selection = "logit")
)

expect_equal(y22_corr_all$output$mean, 0.62883878, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_all$output$SE, 0.069170412, tolerance = 0.0001) ## true value for this sim
expect_true(y22_corr_all$confidence_interval$lower_bound < mean(Y_22) &
              y22_corr_all$confidence_interval$upper_bound > mean(Y_22)) ## conf int

##### all target variables  ---------------------------------------------------------------

expect_silent(
  y_all_corr_all <- nonprob(selection = X_formula,
                            target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                            data = sample_B1,
                            pop_totals = X_totals,
                            method_selection = "logit")
)

expect_identical(y_all_corr_all$output$mean,
                 c(y11_corr_all$output$mean, y12_corr_all$output$mean, y21_corr_all$output$mean, y22_corr_all$output$mean))

expect_identical(y_all_corr_all$output$SE,
                 c(y11_corr_all$output$SE, y12_corr_all$output$SE, y21_corr_all$output$SE, y22_corr_all$output$SE))

expect_identical(y_all_corr_all$confidence_interval,
                 data.frame(lower_bound = c(y11_corr_all$confidence_interval$lower_bound,
                                            y12_corr_all$confidence_interval$lower_bound,
                                            y21_corr_all$confidence_interval$lower_bound,
                                            y22_corr_all$confidence_interval$lower_bound),
                            upper_bound = c(y11_corr_all$confidence_interval$upper_bound,
                                            y12_corr_all$confidence_interval$upper_bound,
                                            y21_corr_all$confidence_interval$upper_bound,
                                            y22_corr_all$confidence_interval$upper_bound),
                            row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))


#### variable selection  ------------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## y_11
expect_silent(
  y11_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_11,
                           data = sample_B1,
                           pop_totals = X_totals[1:11],
                           method_selection = "logit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y11_corr_scad$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_scad$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
expect_false(y11_corr_scad$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_scad$confidence_interval$upper_bound > mean(Y_11)) ## conf int
expect_true(NROW(y11_corr_scad$selection$coefficients) == 2)

## y_12
expect_silent(
  y12_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_12,
                           data = sample_B1,
                           pop_totals = X_totals[1:11],
                           method_selection = "logit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y12_corr_scad$output$mean, 6.9530644, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_scad$output$SE, 0.15341599, tolerance = 0.0001) ## true value for this sim
expect_true(y12_corr_scad$confidence_interval$lower_bound < mean(Y_12) &
               y12_corr_scad$confidence_interval$upper_bound > mean(Y_12)) ## conf int
expect_true(NROW(y12_corr_scad$selection$coefficients) == 2)

## y_21
expect_silent(
  y21_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_21,
                           data = sample_B1,
                           pop_totals = X_totals[1:11],
                           method_selection = "logit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y21_corr_scad$output$mean, 0.78264707, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_scad$output$SE, 0.0090012565, tolerance = 0.0001) ## true value for this sim
expect_false(y21_corr_scad$confidence_interval$lower_bound < mean(Y_21) &
               y21_corr_scad$confidence_interval$upper_bound > mean(Y_21)) ## conf int
expect_true(NROW(y21_corr_scad$selection$coefficients) == 2)

## y_22
expect_silent(
  y22_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_22,
                           data = sample_B1,
                           pop_totals = X_totals[1:11],
                           method_selection = "logit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y22_corr_scad$output$mean, 0.57680653, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_scad$output$SE, 0.011240221, tolerance = 0.0001) ## true value for this sim
expect_false(y22_corr_scad$confidence_interval$lower_bound < mean(Y_22) &
               y22_corr_scad$confidence_interval$upper_bound > mean(Y_22)) ## conf int
expect_true(NROW(y22_corr_scad$selection$coefficients) == 2)


## lasso

# expect_silent(
#   y11_corr_lasso <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_11,
#                            data = sample_B1,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "logit",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "lasso"))
# )
#
# expect_equal(y11_corr_lasso$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_lasso$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
# expect_false(y11_corr_lasso$confidence_interval$lower_bound < mean(Y_11) &
#                y11_corr_lasso$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_lasso$selection$coefficients) == 2)

## MCP

# expect_silent(
#   y11_corr_mcp <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                             target = ~ Y_11,
#                             data = sample_B1,
#                             pop_totals = X_totals[1:11],
#                             method_selection = "logit",
#                             control_inference = controlInf(vars_selection = TRUE),
#                             control_selection = controlSel(penalty = "MCP"))
# )
#
# expect_equal(y11_corr_lasso$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_lasso$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
# expect_false(y11_corr_lasso$confidence_interval$lower_bound < mean(Y_11) &
#                y11_corr_lasso$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_lasso$selection$coefficients) == 2)

##### all target variables  ---------------------------------------------------------------

# expect_silent(
#   y_all_corr_all <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                             target = ~ Y_11, + Y_12 + Y_21 + Y_22,
#                             data = sample_B1,
#                             pop_totals = X_totals[1:11],
#                             method_selection = "logit",
#                             control_inference = controlInf(vars_selection = TRUE),
#                             control_selection = controlSel(penalty = "SCAD", nfolds = 5),
#                             verbose = T)
# )


## non-linear case ------------------------------------------------------------------------
#### correctly specified variables --------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
expect_silent(
  y11_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_11,
                          data = sample_B2,
                          pop_totals = X_totals[1:5],
                          method_selection = "logit")
)

expect_equal(y11_corr_one$output$mean, 1.9842061, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_one$output$SE, 0.1141722, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_one$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_one$confidence_interval$upper_bound > mean(Y_11)) ## conf int


## for y12
expect_silent(
  y12_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_12,
                          data = sample_B2,
                          pop_totals = X_totals[1:5],
                          method_selection = "logit")
)
expect_equal(y12_corr_one$output$mean, 5.7852985, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_one$output$SE, 0.25978825, tolerance = 0.0001) ## true value for this sim
expect_false(y12_corr_one$confidence_interval$lower_bound < mean(Y_12) &
              y12_corr_one$confidence_interval$upper_bound > mean(Y_12)) ## conf int


## for y21
expect_silent(
  y21_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_21,
                          data = sample_B2,
                          pop_totals = X_totals[1:5],
                          method_selection = "logit")
)

expect_equal(y21_corr_one$output$mean, 0.61846098, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_one$output$SE, 0.026561368, tolerance = 0.0001) ## true value for this sim
expect_true(y21_corr_one$confidence_interval$lower_bound < mean(Y_21) &
              y21_corr_one$confidence_interval$upper_bound > mean(Y_21)) ## conf int

## for y22
expect_silent(
  y22_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_22,
                          data = sample_B2,
                          pop_totals = X_totals[1:5],
                          method_selection = "logit")
)

expect_equal(y22_corr_one$output$mean, 0.64512141, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_one$output$SE, 0.026011784, tolerance = 0.0001) ## true value for this sim
expect_true(y22_corr_one$confidence_interval$lower_bound < mean(Y_22) &
              y22_corr_one$confidence_interval$upper_bound > mean(Y_22)) ## conf int

##### all target variables  ---------------------------------------------------------------

expect_silent(
  y_all_corr <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                        target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                        data = sample_B2,
                        pop_totals = X_totals[1:5],
                        method_selection = "logit")
)

expect_identical(y_all_corr$output$mean,
                 c(y11_corr_one$output$mean, y12_corr_one$output$mean,
                   y21_corr_one$output$mean, y22_corr_one$output$mean))

expect_identical(y_all_corr$output$SE,
                 c(y11_corr_one$output$SE, y12_corr_one$output$SE,
                   y21_corr_one$output$SE, y22_corr_one$output$SE))

expect_identical(y_all_corr$confidence_interval,
                 data.frame(lower_bound = c(y11_corr_one$confidence_interval$lower_bound,
                                            y12_corr_one$confidence_interval$lower_bound,
                                            y21_corr_one$confidence_interval$lower_bound,
                                            y22_corr_one$confidence_interval$lower_bound),
                            upper_bound = c(y11_corr_one$confidence_interval$upper_bound,
                                            y12_corr_one$confidence_interval$upper_bound,
                                            y21_corr_one$confidence_interval$upper_bound,
                                            y22_corr_one$confidence_interval$upper_bound),
                            row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))





#### all X variables variables ------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
expect_silent(
  y11_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_11,
                          data = sample_B2,
                          pop_totals = X_totals,
                          method_selection = "logit")
)

expect_equal(y11_corr_all$output$mean, 1.9820754, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_all$output$SE, 0.13574794, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_all$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_all$confidence_interval$upper_bound > mean(Y_11)) ## conf int


## for y12
expect_silent(
  y12_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_12,
                          data = sample_B2,
                          pop_totals = X_totals,
                          method_selection = "logit")
)
expect_equal(y12_corr_all$output$mean, 5.6417776, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_all$output$SE, 0.27204781, tolerance = 0.0001) ## true value for this sim
expect_false(y12_corr_all$confidence_interval$lower_bound < mean(Y_12) &
              y12_corr_all$confidence_interval$upper_bound > mean(Y_12)) ## conf int


## for y21
expect_silent(
  y21_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_21,
                          data = sample_B2,
                          pop_totals = X_totals,
                          method_selection = "logit")
)

expect_equal(y21_corr_all$output$mean, 0.61556545, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_all$output$SE, 0.028680898, tolerance = 0.0001) ## true value for this sim
expect_true(y21_corr_all$confidence_interval$lower_bound < mean(Y_21) &
              y21_corr_all$confidence_interval$upper_bound > mean(Y_21)) ## conf int

## for y22
expect_silent(
  y22_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_22,
                          data = sample_B2,
                          pop_totals = X_totals,
                          method_selection = "logit")
)

expect_equal(y22_corr_all$output$mean, 0.57173102, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_all$output$SE, 0.024335536, tolerance = 0.0001) ## true value for this sim
expect_false(y22_corr_all$confidence_interval$lower_bound < mean(Y_22) &
              y22_corr_all$confidence_interval$upper_bound > mean(Y_22)) ## conf int

##### all target variables  ---------------------------------------------------------------

expect_silent(
  y_all_corr_all <- nonprob(selection = X_formula,
                            target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                            data = sample_B2,
                            pop_totals = X_totals,
                            method_selection = "logit")
)

expect_identical(y_all_corr_all$output$mean,
                 c(y11_corr_all$output$mean, y12_corr_all$output$mean, y21_corr_all$output$mean, y22_corr_all$output$mean))

expect_identical(y_all_corr_all$output$SE,
                 c(y11_corr_all$output$SE, y12_corr_all$output$SE, y21_corr_all$output$SE, y22_corr_all$output$SE))

expect_identical(y_all_corr_all$confidence_interval,
                 data.frame(lower_bound = c(y11_corr_all$confidence_interval$lower_bound,
                                            y12_corr_all$confidence_interval$lower_bound,
                                            y21_corr_all$confidence_interval$lower_bound,
                                            y22_corr_all$confidence_interval$lower_bound),
                            upper_bound = c(y11_corr_all$confidence_interval$upper_bound,
                                            y12_corr_all$confidence_interval$upper_bound,
                                            y21_corr_all$confidence_interval$upper_bound,
                                            y22_corr_all$confidence_interval$upper_bound),
                            row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))



#### variable selection  ------------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## y_11
expect_silent(
  y11_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_11,
                           data = sample_B2,
                           pop_totals = X_totals[1:11],
                           method_selection = "logit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y11_corr_scad$output$mean, 1.8810431, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_scad$output$SE, 0.059381217, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_scad$confidence_interval$lower_bound < mean(Y_11) &
               y11_corr_scad$confidence_interval$upper_bound > mean(Y_11)) ## conf int
expect_true(NROW(y11_corr_scad$selection$coefficients) == 2)

## y_12
expect_silent(
  y12_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_12,
                           data = sample_B2,
                           pop_totals = X_totals[1:11],
                           method_selection = "logit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y12_corr_scad$output$mean, 5.796713, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_scad$output$SE, 0.14583111, tolerance = 0.0001) ## true value for this sim
expect_false(y12_corr_scad$confidence_interval$lower_bound < mean(Y_12) &
               y12_corr_scad$confidence_interval$upper_bound > mean(Y_12)) ## conf int
expect_true(NROW(y12_corr_scad$selection$coefficients) == 2)

## y_21
expect_silent(
  y21_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_21,
                           data = sample_B2,
                           pop_totals = X_totals[1:11],
                           method_selection = "logit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y21_corr_scad$output$mean, 0.6060074, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_scad$output$SE, 0.010194911, tolerance = 0.0001) ## true value for this sim
expect_false(y21_corr_scad$confidence_interval$lower_bound < mean(Y_21) &
               y21_corr_scad$confidence_interval$upper_bound > mean(Y_21)) ## conf int
expect_true(NROW(y21_corr_scad$selection$coefficients) == 2)

## y_22
expect_silent(
  y22_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_22,
                           data = sample_B2,
                           pop_totals = X_totals[1:11],
                           method_selection = "logit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y22_corr_scad$output$mean, 0.64707641, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_scad$output$SE, 0.0099648982, tolerance = 0.0001) ## true value for this sim
expect_true(y22_corr_scad$confidence_interval$lower_bound < mean(Y_22) &
               y22_corr_scad$confidence_interval$upper_bound > mean(Y_22)) ## conf int
expect_true(NROW(y22_corr_scad$selection$coefficients) == 2)


## lasso

# expect_silent(
#   y11_corr_lasso <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_11,
#                            data = sample_B1,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "logit",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "lasso"))
# )
#
# expect_equal(y11_corr_lasso$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_lasso$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
# expect_false(y11_corr_lasso$confidence_interval$lower_bound < mean(Y_11) &
#                y11_corr_lasso$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_lasso$selection$coefficients) == 2)

## MCP

# expect_silent(
#   y11_corr_mcp <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                             target = ~ Y_11,
#                             data = sample_B1,
#                             pop_totals = X_totals[1:11],
#                             method_selection = "logit",
#                             control_inference = controlInf(vars_selection = TRUE),
#                             control_selection = controlSel(penalty = "MCP"))
# )
#
# expect_equal(y11_corr_lasso$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_lasso$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
# expect_false(y11_corr_lasso$confidence_interval$lower_bound < mean(Y_11) &
#                y11_corr_lasso$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_lasso$selection$coefficients) == 2)

##### all target variables  ---------------------------------------------------------------

# expect_silent(
#   y_all_corr_all <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                             target = ~ Y_11, + Y_12 + Y_21 + Y_22,
#                             data = sample_B1,
#                             pop_totals = X_totals[1:11],
#                             method_selection = "logit",
#                             control_inference = controlInf(vars_selection = TRUE),
#                             control_selection = controlSel(penalty = "SCAD", nfolds = 5),
#                             verbose = T)
# )




# check probit ----------------------------------------------------------------------------
## linear case ----------------------------------------------------------------------------
#### correctly specified variables --------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
expect_silent(
  y11_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_11,
                          data = sample_B1,
                          pop_totals = X_totals[1:5],
                          method_selection = "probit")
)

expect_equal(y11_corr_one$output$mean, 2.1985815, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_one$output$SE, 0.12852317, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_one$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_one$confidence_interval$upper_bound > mean(Y_11)) ## conf int


## for y12
expect_silent(
  y12_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_12,
                          data = sample_B1,
                          pop_totals = X_totals[1:5],
                          method_selection = "probit")
)
expect_equal(y12_corr_one$output$mean, 7.1226154, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_one$output$SE, 0.83353647, tolerance = 0.0001) ## true value for this sim
expect_true(y12_corr_one$confidence_interval$lower_bound < mean(Y_12) &
              y12_corr_one$confidence_interval$upper_bound > mean(Y_12)) ## conf int


## for y21
expect_silent(
  y21_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_21,
                          data = sample_B1,
                          pop_totals = X_totals[1:5],
                          method_selection = "probit")
)

expect_equal(y21_corr_one$output$mean, 0.65778764, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_one$output$SE, 0.066138969, tolerance = 0.0001) ## true value for this sim
expect_true(y21_corr_one$confidence_interval$lower_bound < mean(Y_21) &
              y21_corr_one$confidence_interval$upper_bound > mean(Y_21)) ## conf int

## for y22
expect_silent(
  y22_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_22,
                          data = sample_B1,
                          pop_totals = X_totals[1:5],
                          method_selection = "probit")
)

expect_equal(y22_corr_one$output$mean, 0.58148029, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_one$output$SE, 0.068800601, tolerance = 0.0001) ## true value for this sim
expect_true(y22_corr_one$confidence_interval$lower_bound < mean(Y_22) &
              y22_corr_one$confidence_interval$upper_bound > mean(Y_22)) ## conf int

##### all target variables  ---------------------------------------------------------------

expect_silent(
  y_all_corr <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                        target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                        data = sample_B1,
                        pop_totals = X_totals[1:5],
                        method_selection = "probit")
)

expect_identical(y_all_corr$output$mean,
                 c(y11_corr_one$output$mean, y12_corr_one$output$mean,
                   y21_corr_one$output$mean, y22_corr_one$output$mean))

expect_identical(y_all_corr$output$SE,
                 c(y11_corr_one$output$SE, y12_corr_one$output$SE,
                   y21_corr_one$output$SE, y22_corr_one$output$SE))

expect_identical(y_all_corr$confidence_interval,
                 data.frame(lower_bound = c(y11_corr_one$confidence_interval$lower_bound,
                                            y12_corr_one$confidence_interval$lower_bound,
                                            y21_corr_one$confidence_interval$lower_bound,
                                            y22_corr_one$confidence_interval$lower_bound),
                            upper_bound = c(y11_corr_one$confidence_interval$upper_bound,
                                            y12_corr_one$confidence_interval$upper_bound,
                                            y21_corr_one$confidence_interval$upper_bound,
                                            y22_corr_one$confidence_interval$upper_bound),
                            row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))





#### all X variables variables ------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
expect_silent(
  y11_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_11,
                          data = sample_B1,
                          pop_totals = X_totals,
                          method_selection = "probit")
)

expect_equal(y11_corr_all$output$mean, 2.0102108, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_all$output$SE, 0.16422236, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_all$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_all$confidence_interval$upper_bound > mean(Y_11)) ## conf int


## for y12
expect_silent(
  y12_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_12,
                          data = sample_B1,
                          pop_totals = X_totals,
                          method_selection = "probit")
)
expect_equal(y12_corr_all$output$mean, 6.6878138, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_all$output$SE, 0.60740161, tolerance = 0.0001) ## true value for this sim
expect_true(y12_corr_all$confidence_interval$lower_bound < mean(Y_12) &
              y12_corr_all$confidence_interval$upper_bound > mean(Y_12)) ## conf int


## for y21
expect_silent(
  y21_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_21,
                          data = sample_B1,
                          pop_totals = X_totals,
                          method_selection = "probit")
)

expect_equal(y21_corr_all$output$mean, 0.64617883, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_all$output$SE, 0.074037561, tolerance = 0.0001) ## true value for this sim
expect_true(y21_corr_all$confidence_interval$lower_bound < mean(Y_21) &
              y21_corr_all$confidence_interval$upper_bound > mean(Y_21)) ## conf int

## for y22
expect_silent(
  y22_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_22,
                          data = sample_B1,
                          pop_totals = X_totals,
                          method_selection = "probit")
)

expect_equal(y22_corr_all$output$mean, 0.62456575, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_all$output$SE, 0.083014906, tolerance = 0.0001) ## true value for this sim
expect_true(y22_corr_all$confidence_interval$lower_bound < mean(Y_22) &
              y22_corr_all$confidence_interval$upper_bound > mean(Y_22)) ## conf int

##### all target variables  ---------------------------------------------------------------

expect_silent(
  y_all_corr_all <- nonprob(selection = X_formula,
                            target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                            data = sample_B1,
                            pop_totals = X_totals,
                            method_selection = "probit")
)

expect_identical(y_all_corr_all$output$mean,
                 c(y11_corr_all$output$mean, y12_corr_all$output$mean, y21_corr_all$output$mean, y22_corr_all$output$mean))

expect_identical(y_all_corr_all$output$SE,
                 c(y11_corr_all$output$SE, y12_corr_all$output$SE, y21_corr_all$output$SE, y22_corr_all$output$SE))

expect_identical(y_all_corr_all$confidence_interval,
                 data.frame(lower_bound = c(y11_corr_all$confidence_interval$lower_bound,
                                            y12_corr_all$confidence_interval$lower_bound,
                                            y21_corr_all$confidence_interval$lower_bound,
                                            y22_corr_all$confidence_interval$lower_bound),
                            upper_bound = c(y11_corr_all$confidence_interval$upper_bound,
                                            y12_corr_all$confidence_interval$upper_bound,
                                            y21_corr_all$confidence_interval$upper_bound,
                                            y22_corr_all$confidence_interval$upper_bound),
                            row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))



#### variable selection  ------------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## y_11
expect_silent(
  y11_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_11,
                           data = sample_B1,
                           pop_totals = X_totals[1:11],
                           method_selection = "probit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y11_corr_scad$output$mean, 3.0633399, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_scad$output$SE, 0.049384577, tolerance = 0.0001) ## true value for this sim
expect_false(y11_corr_scad$confidence_interval$lower_bound < mean(Y_11) &
               y11_corr_scad$confidence_interval$upper_bound > mean(Y_11)) ## conf int
expect_true(NROW(y11_corr_scad$selection$coefficients) == 2)

## y_12
expect_silent(
  y12_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_12,
                           data = sample_B1,
                           pop_totals = X_totals[1:11],
                           method_selection = "probit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y12_corr_scad$output$mean, 6.9420676, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_scad$output$SE, 0.15578019, tolerance = 0.0001) ## true value for this sim
expect_true(y12_corr_scad$confidence_interval$lower_bound < mean(Y_12) &
               y12_corr_scad$confidence_interval$upper_bound > mean(Y_12)) ## conf int
expect_true(NROW(y12_corr_scad$selection$coefficients) == 2)

## y_21
expect_silent(
  y21_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_21,
                           data = sample_B1,
                           pop_totals = X_totals[1:11],
                           method_selection = "probit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y21_corr_scad$output$mean, 0.78324543, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_scad$output$SE, 0.0090654702, tolerance = 0.0001) ## true value for this sim
expect_false(y21_corr_scad$confidence_interval$lower_bound < mean(Y_21) &
               y21_corr_scad$confidence_interval$upper_bound > mean(Y_21)) ## conf int
expect_true(NROW(y21_corr_scad$selection$coefficients) == 2)

## y_22
expect_silent(
  y22_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_22,
                           data = sample_B1,
                           pop_totals = X_totals[1:11],
                           method_selection = "probit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y22_corr_scad$output$mean, 0.57672297, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_scad$output$SE, 0.011433876, tolerance = 0.0001) ## true value for this sim
expect_false(y22_corr_scad$confidence_interval$lower_bound < mean(Y_22) &
               y22_corr_scad$confidence_interval$upper_bound > mean(Y_22)) ## conf int
expect_true(NROW(y22_corr_scad$selection$coefficients) == 2)


## lasso

# expect_silent(
#   y11_corr_lasso <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_11,
#                            data = sample_B1,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "probit",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "lasso"))
# )
#
# expect_equal(y11_corr_lasso$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_lasso$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
# expect_false(y11_corr_lasso$confidence_interval$lower_bound < mean(Y_11) &
#                y11_corr_lasso$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_lasso$selection$coefficients) == 2)

## MCP

# expect_silent(
#   y11_corr_mcp <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                             target = ~ Y_11,
#                             data = sample_B1,
#                             pop_totals = X_totals[1:11],
#                             method_selection = "probit",
#                             control_inference = controlInf(vars_selection = TRUE),
#                             control_selection = controlSel(penalty = "MCP"))
# )
#
# expect_equal(y11_corr_lasso$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_lasso$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
# expect_false(y11_corr_lasso$confidence_interval$lower_bound < mean(Y_11) &
#                y11_corr_lasso$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_lasso$selection$coefficients) == 2)

##### all target variables  ---------------------------------------------------------------

# expect_silent(
#   y_all_corr_all <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                             target = ~ Y_11, + Y_12 + Y_21 + Y_22,
#                             data = sample_B1,
#                             pop_totals = X_totals[1:11],
#                             method_selection = "probit",
#                             control_inference = controlInf(vars_selection = TRUE),
#                             control_selection = controlSel(penalty = "SCAD", nfolds = 5),
#                             verbose = T)
# )



## non-linear case ------------------------------------------------------------------------
#### correctly specified variables --------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
expect_silent(
  y11_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_11,
                          data = sample_B2,
                          pop_totals = X_totals[1:5],
                          method_selection = "probit")
)

expect_equal(y11_corr_one$output$mean, 1.984232, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_one$output$SE, 0.10478658, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_one$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_one$confidence_interval$upper_bound > mean(Y_11)) ## conf int


## for y12
expect_silent(
  y12_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_12,
                          data = sample_B2,
                          pop_totals = X_totals[1:5],
                          method_selection = "probit")
)
expect_equal(y12_corr_one$output$mean, 5.7851563, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_one$output$SE, 0.20174126, tolerance = 0.0001) ## true value for this sim
expect_false(y12_corr_one$confidence_interval$lower_bound < mean(Y_12) &
               y12_corr_one$confidence_interval$upper_bound > mean(Y_12)) ## conf int


## for y21
expect_silent(
  y21_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_21,
                          data = sample_B2,
                          pop_totals = X_totals[1:5],
                          method_selection = "probit")
)

expect_equal(y21_corr_one$output$mean, 0.61846032, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_one$output$SE, 0.031103296, tolerance = 0.0001) ## true value for this sim
expect_true(y21_corr_one$confidence_interval$lower_bound < mean(Y_21) &
              y21_corr_one$confidence_interval$upper_bound > mean(Y_21)) ## conf int

## for y22
expect_silent(
  y22_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                          target = ~ Y_22,
                          data = sample_B2,
                          pop_totals = X_totals[1:5],
                          method_selection = "probit")
)

expect_equal(y22_corr_one$output$mean, 0.64509743, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_one$output$SE, 0.030385359, tolerance = 0.0001) ## true value for this sim
expect_true(y22_corr_one$confidence_interval$lower_bound < mean(Y_22) &
              y22_corr_one$confidence_interval$upper_bound > mean(Y_22)) ## conf int

##### all target variables  ---------------------------------------------------------------

expect_silent(
  y_all_corr <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                        target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                        data = sample_B2,
                        pop_totals = X_totals[1:5],
                        method_selection = "probit")
)

expect_identical(y_all_corr$output$mean,
                 c(y11_corr_one$output$mean, y12_corr_one$output$mean,
                   y21_corr_one$output$mean, y22_corr_one$output$mean))

expect_identical(y_all_corr$output$SE,
                 c(y11_corr_one$output$SE, y12_corr_one$output$SE,
                   y21_corr_one$output$SE, y22_corr_one$output$SE))

expect_identical(y_all_corr$confidence_interval,
                 data.frame(lower_bound = c(y11_corr_one$confidence_interval$lower_bound,
                                            y12_corr_one$confidence_interval$lower_bound,
                                            y21_corr_one$confidence_interval$lower_bound,
                                            y22_corr_one$confidence_interval$lower_bound),
                            upper_bound = c(y11_corr_one$confidence_interval$upper_bound,
                                            y12_corr_one$confidence_interval$upper_bound,
                                            y21_corr_one$confidence_interval$upper_bound,
                                            y22_corr_one$confidence_interval$upper_bound),
                            row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))





#### all X variables variables ------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
expect_silent(
  y11_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_11,
                          data = sample_B2,
                          pop_totals = X_totals,
                          method_selection = "probit")
)

expect_equal(y11_corr_all$output$mean, 1.981833, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_all$output$SE, 0.12730647, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_all$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_all$confidence_interval$upper_bound > mean(Y_11)) ## conf int


## for y12
expect_silent(
  y12_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_12,
                          data = sample_B2,
                          pop_totals = X_totals,
                          method_selection = "probit")
)
expect_equal(y12_corr_all$output$mean, 5.6436133, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_all$output$SE, 0.21772187, tolerance = 0.0001) ## true value for this sim
expect_false(y12_corr_all$confidence_interval$lower_bound < mean(Y_12) &
               y12_corr_all$confidence_interval$upper_bound > mean(Y_12)) ## conf int


## for y21
expect_silent(
  y21_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_21,
                          data = sample_B2,
                          pop_totals = X_totals,
                          method_selection = "probit")
)

expect_equal(y21_corr_all$output$mean, 0.61562293, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_all$output$SE, 0.033259168, tolerance = 0.0001) ## true value for this sim
expect_true(y21_corr_all$confidence_interval$lower_bound < mean(Y_21) &
              y21_corr_all$confidence_interval$upper_bound > mean(Y_21)) ## conf int

## for y22
expect_silent(
  y22_corr_all <- nonprob(selection = X_formula,
                          target = ~ Y_22,
                          data = sample_B2,
                          pop_totals = X_totals,
                          method_selection = "probit")
)

expect_equal(y22_corr_all$output$mean, 0.57173464, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_all$output$SE, 0.028962455, tolerance = 0.0001) ## true value for this sim
expect_false(y22_corr_all$confidence_interval$lower_bound < mean(Y_22) &
               y22_corr_all$confidence_interval$upper_bound > mean(Y_22)) ## conf int

##### all target variables  ---------------------------------------------------------------

expect_silent(
  y_all_corr_all <- nonprob(selection = X_formula,
                            target = ~ Y_11 + Y_12 + Y_21 + Y_22,
                            data = sample_B2,
                            pop_totals = X_totals,
                            method_selection = "probit")
)

expect_identical(y_all_corr_all$output$mean,
                 c(y11_corr_all$output$mean, y12_corr_all$output$mean,
                   y21_corr_all$output$mean, y22_corr_all$output$mean))

expect_identical(y_all_corr_all$output$SE,
                 c(y11_corr_all$output$SE, y12_corr_all$output$SE,
                   y21_corr_all$output$SE, y22_corr_all$output$SE))

expect_identical(y_all_corr_all$confidence_interval,
                 data.frame(lower_bound = c(y11_corr_all$confidence_interval$lower_bound,
                                            y12_corr_all$confidence_interval$lower_bound,
                                            y21_corr_all$confidence_interval$lower_bound,
                                            y22_corr_all$confidence_interval$lower_bound),
                            upper_bound = c(y11_corr_all$confidence_interval$upper_bound,
                                            y12_corr_all$confidence_interval$upper_bound,
                                            y21_corr_all$confidence_interval$upper_bound,
                                            y22_corr_all$confidence_interval$upper_bound),
                            row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))



#### variable selection  ------------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## y_11
expect_silent(
  y11_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_11,
                           data = sample_B2,
                           pop_totals = X_totals[1:11],
                           method_selection = "probit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y11_corr_scad$output$mean, 1.8810431, tolerance = 0.0001) ## true value for this sim
expect_equal(y11_corr_scad$output$SE, 0.059381198, tolerance = 0.0001) ## true value for this sim
expect_true(y11_corr_scad$confidence_interval$lower_bound < mean(Y_11) &
              y11_corr_scad$confidence_interval$upper_bound > mean(Y_11)) ## conf int
expect_true(NROW(y11_corr_scad$selection$coefficients) == 2)

## y_12
expect_silent(
  y12_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_12,
                           data = sample_B2,
                           pop_totals = X_totals[1:11],
                           method_selection = "probit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y12_corr_scad$output$mean, 5.7967136, tolerance = 0.0001) ## true value for this sim
expect_equal(y12_corr_scad$output$SE, 0.14583128, tolerance = 0.0001) ## true value for this sim
expect_false(y12_corr_scad$confidence_interval$lower_bound < mean(Y_12) &
               y12_corr_scad$confidence_interval$upper_bound > mean(Y_12)) ## conf int
expect_true(NROW(y12_corr_scad$selection$coefficients) == 2)

## y_21
expect_silent(
  y21_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_21,
                           data = sample_B2,
                           pop_totals = X_totals[1:11],
                           method_selection = "probit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y21_corr_scad$output$mean, 0.60600756, tolerance = 0.0001) ## true value for this sim
expect_equal(y21_corr_scad$output$SE, 0.010194928, tolerance = 0.0001) ## true value for this sim
expect_false(y21_corr_scad$confidence_interval$lower_bound < mean(Y_21) &
               y21_corr_scad$confidence_interval$upper_bound > mean(Y_21)) ## conf int
expect_true(NROW(y21_corr_scad$selection$coefficients) == 2)

## y_22
expect_silent(
  y22_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                           target = ~ Y_22,
                           data = sample_B2,
                           pop_totals = X_totals[1:11],
                           method_selection = "probit",
                           control_inference = controlInf(vars_selection = TRUE),
                           control_selection = controlSel(penalty = "SCAD", nfolds = 5))
)

expect_equal(y22_corr_scad$output$mean, 0.64707626, tolerance = 0.0001) ## true value for this sim
expect_equal(y22_corr_scad$output$SE, 0.0099648984, tolerance = 0.0001) ## true value for this sim
expect_true(y22_corr_scad$confidence_interval$lower_bound < mean(Y_22) &
               y22_corr_scad$confidence_interval$upper_bound > mean(Y_22)) ## conf int
expect_true(NROW(y22_corr_scad$selection$coefficients) == 2)


## lasso

# expect_silent(
#   y11_corr_lasso <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_11,
#                            data = sample_B1,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "probit",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "lasso"))
# )
#
# expect_equal(y11_corr_lasso$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_lasso$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
# expect_false(y11_corr_lasso$confidence_interval$lower_bound < mean(Y_11) &
#                y11_corr_lasso$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_lasso$selection$coefficients) == 2)

## MCP

# expect_silent(
#   y11_corr_mcp <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                             target = ~ Y_11,
#                             data = sample_B1,
#                             pop_totals = X_totals[1:11],
#                             method_selection = "probit",
#                             control_inference = controlInf(vars_selection = TRUE),
#                             control_selection = controlSel(penalty = "MCP"))
# )
#
# expect_equal(y11_corr_lasso$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_lasso$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
# expect_false(y11_corr_lasso$confidence_interval$lower_bound < mean(Y_11) &
#                y11_corr_lasso$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_lasso$selection$coefficients) == 2)

##### all target variables  ---------------------------------------------------------------

# expect_silent(
#   y_all_corr_all <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                             target = ~ Y_11, + Y_12 + Y_21 + Y_22,
#                             data = sample_B1,
#                             pop_totals = X_totals[1:11],
#                             method_selection = "probit",
#                             control_inference = controlInf(vars_selection = TRUE),
#                             control_selection = controlSel(penalty = "SCAD", nfolds = 5),
#                             verbose = T)
# )





# check cloglog -----------------------------------------------------
## linear case ----------------------------------------------------------------------------
#### correctly specified variables --------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
# expect_silent(
#   y11_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                           target = ~ Y_11,
#                           data = sample_B1,
#                           pop_totals = X_totals[1:5],
#                           method_selection = "cloglog")
# )
#
# expect_equal(y11_corr_one$output$mean, 2.015028, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_one$output$SE, 0.1713938, tolerance = 0.0001) ## true value for this sim
# expect_true(y11_corr_one$confidence_interval$lower_bound < mean(Y_11) &
#               y11_corr_one$confidence_interval$upper_bound > mean(Y_11)) ## conf int
#
#
# ## for y12
# expect_silent(
#   y12_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                           target = ~ Y_12,
#                           data = sample_B1,
#                           pop_totals = X_totals[1:5],
#                           method_selection = "cloglog")
# )
# expect_equal(y12_corr_one$output$mean, 6.617158, tolerance = 0.0001) ## true value for this sim
# expect_equal(y12_corr_one$output$SE, 0.6183295, tolerance = 0.0001) ## true value for this sim
# expect_true(y12_corr_one$confidence_interval$lower_bound < mean(Y_12) &
#               y12_corr_one$confidence_interval$upper_bound > mean(Y_12)) ## conf int
#
#
# ## for y21
# expect_silent(
#   y21_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                           target = ~ Y_21,
#                           data = sample_B1,
#                           pop_totals = X_totals[1:5],
#                           method_selection = "cloglog")
# )
#
# expect_equal(y21_corr_one$output$mean, 0.6589795, tolerance = 0.0001) ## true value for this sim
# expect_equal(y21_corr_one$output$SE, 0.05471915, tolerance = 0.0001) ## true value for this sim
# expect_true(y21_corr_one$confidence_interval$lower_bound < mean(Y_21) &
#               y21_corr_one$confidence_interval$upper_bound > mean(Y_21)) ## conf int
#
# ## for y22
# expect_silent(
#   y22_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                           target = ~ Y_22,
#                           data = sample_B1,
#                           pop_totals = X_totals[1:5],
#                           method_selection = "cloglog")
# )
#
# expect_equal(y22_corr_one$output$mean, 0.6872036, tolerance = 0.0001) ## true value for this sim
# expect_equal(y22_corr_one$output$SE, 0.06540456, tolerance = 0.0001) ## true value for this sim
# expect_true(y22_corr_one$confidence_interval$lower_bound < mean(Y_22) &
#               y22_corr_one$confidence_interval$upper_bound > mean(Y_22)) ## conf int
#
# ##### all target variables  ---------------------------------------------------------------
#
# expect_silent(
#   y_all_corr <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                         target = ~ Y_11 + Y_12 + Y_21 + Y_22,
#                         data = sample_B1,
#                         pop_totals = X_totals[1:5],
#                         method_selection = "cloglog")
# )
#
# expect_identical(y_all_corr$output$mean,
#                  c(y11_corr_one$output$mean, y12_corr_one$output$mean,
#                    y21_corr_one$output$mean, y22_corr_one$output$mean))
#
# expect_identical(y_all_corr$output$SE,
#                  c(y11_corr_one$output$SE, y12_corr_one$output$SE,
#                    y21_corr_one$output$SE, y22_corr_one$output$SE))
# #
# expect_identical(y_all_corr$confidence_interval,
#                  data.frame(lower_bound = c(y11_corr_one$confidence_interval$lower_bound,
#                                             y12_corr_one$confidence_interval$lower_bound,
#                                             y21_corr_one$confidence_interval$lower_bound,
#                                             y22_corr_one$confidence_interval$lower_bound),
#                             upper_bound = c(y11_corr_one$confidence_interval$upper_bound,
#                                             y12_corr_one$confidence_interval$upper_bound,
#                                             y21_corr_one$confidence_interval$upper_bound,
#                                             y22_corr_one$confidence_interval$upper_bound),
#                             row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))
#
#
#
# #### all X variables variables ------------------------------------------------------------
# ##### one target variable  ----------------------------------------------------------------
#
# ## for y11
# expect_silent(
#   y11_corr_all <- nonprob(selection = X_formula,
#                           target = ~ Y_11,
#                           data = sample_B1,
#                           pop_totals = X_totals,
#                           method_selection = "cloglog")
# )
#
# expect_equal(y11_corr_all$output$mean, 1.956418, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_all$output$SE, 0.2019104, tolerance = 0.0001) ## true value for this sim
# expect_true(y11_corr_all$confidence_interval$lower_bound < mean(Y_11) &
#               y11_corr_all$confidence_interval$upper_bound > mean(Y_11)) ## conf int
#
#
# ## for y12
# expect_silent(
#   y12_corr_all <- nonprob(selection = X_formula,
#                           target = ~ Y_12,
#                           data = sample_B1,
#                           pop_totals = X_totals,
#                           method_selection = "cloglog")
# )
# expect_equal(y12_corr_all$output$mean, 6.571494, tolerance = 0.0001) ## true value for this sim
# expect_equal(y12_corr_all$output$SE, 0.7128783, tolerance = 0.0001) ## true value for this sim
# expect_true(y12_corr_all$confidence_interval$lower_bound < mean(Y_12) &
#               y12_corr_all$confidence_interval$upper_bound > mean(Y_12)) ## conf int
#
#
# ## for y21
# expect_silent(
#   y21_corr_all <- nonprob(selection = X_formula,
#                           target = ~ Y_21,
#                           data = sample_B1,
#                           pop_totals = X_totals,
#                           method_selection = "cloglog")
# )
#
# expect_equal(y21_corr_all$output$mean, 0.6496145, tolerance = 0.0001) ## true value for this sim
# expect_equal(y21_corr_all$output$SE, 0.0575112, tolerance = 0.0001) ## true value for this sim
# expect_true(y21_corr_all$confidence_interval$lower_bound < mean(Y_21) &
#               y21_corr_all$confidence_interval$upper_bound > mean(Y_21)) ## conf int
#
# ## for y22
# expect_silent(
#   y22_corr_all <- nonprob(selection = X_formula,
#                           target = ~ Y_22,
#                           data = sample_B1,
#                           pop_totals = X_totals,
#                           method_selection = "cloglog")
# )
#
# expect_equal(y22_corr_all$output$mean, 0.6763156, tolerance = 0.0001) ## true value for this sim
# expect_equal(y22_corr_all$output$SE, 0.06510848, tolerance = 0.0001) ## true value for this sim
# expect_true(y22_corr_all$confidence_interval$lower_bound < mean(Y_22) &
#               y22_corr_all$confidence_interval$upper_bound > mean(Y_22)) ## conf int
#
# ##### all target variables  ---------------------------------------------------------------
#
# expect_silent(
#   y_all_corr_all <- nonprob(selection = X_formula,
#                             target = ~ Y_11 + Y_12 + Y_21 + Y_22,
#                             data = sample_B1,
#                             pop_totals = X_totals,
#                             method_selection = "cloglog")
# )
#
# expect_identical(y_all_corr_all$output$mean,
#                  c(y11_corr_all$output$mean, y12_corr_all$output$mean, y21_corr_all$output$mean, y22_corr_all$output$mean))
#
# expect_identical(y_all_corr_all$output$SE,
#                  c(y11_corr_all$output$SE, y12_corr_all$output$SE, y21_corr_all$output$SE, y22_corr_all$output$SE))
#
# expect_identical(y_all_corr_all$confidence_interval,
#                  data.frame(lower_bound = c(y11_corr_all$confidence_interval$lower_bound,
#                                             y12_corr_all$confidence_interval$lower_bound,
#                                             y21_corr_all$confidence_interval$lower_bound,
#                                             y22_corr_all$confidence_interval$lower_bound),
#                             upper_bound = c(y11_corr_all$confidence_interval$upper_bound,
#                                             y12_corr_all$confidence_interval$upper_bound,
#                                             y21_corr_all$confidence_interval$upper_bound,
#                                             y22_corr_all$confidence_interval$upper_bound),
#                             row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))
#

#### variable selection  ------------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## y_11 to fix
# expect_silent(
#   y11_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_11,
#                            data = sample_B1,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "cloglog",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "SCAD", nfolds = 5))
# )

# expect_equal(y11_corr_scad$output$mean, 3.063926, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_scad$output$SE, 0.04853563, tolerance = 0.0001) ## true value for this sim
# expect_false(y11_corr_scad$confidence_interval$lower_bound < mean(Y_11) &
#                y11_corr_scad$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_scad$selection$coefficients) == 2)

## y_12 to fix
# expect_silent(
#   y12_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_12,
#                            data = sample_B1,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "cloglog",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "SCAD", nfolds = 5))
# )

# expect_equal(y12_corr_scad$output$mean, 6.9530644, tolerance = 0.0001) ## true value for this sim
# expect_equal(y12_corr_scad$output$SE, 0.15341599, tolerance = 0.0001) ## true value for this sim
# expect_true(y12_corr_scad$confidence_interval$lower_bound < mean(Y_12) &
#               y12_corr_scad$confidence_interval$upper_bound > mean(Y_12)) ## conf int
# expect_true(NROW(y12_corr_scad$selection$coefficients) == 2)

## y_21 to fix
# expect_silent(
#   y21_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_21,
#                            data = sample_B1,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "cloglog",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "SCAD", nfolds = 5))
# )

# expect_equal(y21_corr_scad$output$mean, 0.78264707, tolerance = 0.0001) ## true value for this sim
# expect_equal(y21_corr_scad$output$SE, 0.0090012565, tolerance = 0.0001) ## true value for this sim
# expect_false(y21_corr_scad$confidence_interval$lower_bound < mean(Y_21) &
#                y21_corr_scad$confidence_interval$upper_bound > mean(Y_21)) ## conf int
# expect_true(NROW(y21_corr_scad$selection$coefficients) == 2)

## y_22 to fix
# expect_silent(
#   y22_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_22,
#                            data = sample_B1,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "cloglog",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "SCAD", nfolds = 5))
# )

# expect_equal(y22_corr_scad$output$mean, 0.57680653, tolerance = 0.0001) ## true value for this sim
# expect_equal(y22_corr_scad$output$SE, 0.011240221, tolerance = 0.0001) ## true value for this sim
# expect_false(y22_corr_scad$confidence_interval$lower_bound < mean(Y_22) &
#                y22_corr_scad$confidence_interval$upper_bound > mean(Y_22)) ## conf int
# expect_true(NROW(y22_corr_scad$selection$coefficients) == 2)


## non-linear case ------------------------------------------------------------------------
#### correctly specified variables --------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## for y11
# expect_silent(
#   y11_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                           target = ~ Y_11,
#                           data = sample_B2,
#                           pop_totals = X_totals[1:5],
#                           method_selection = "cloglog")
# )
#
# expect_equal(y11_corr_one$output$mean, 2.10861, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_one$output$SE, 0.1243684, tolerance = 0.0001) ## true value for this sim
# expect_true(y11_corr_one$confidence_interval$lower_bound < mean(Y_11) &
#               y11_corr_one$confidence_interval$upper_bound > mean(Y_11)) ## conf int
#
#
# ## for y12
# expect_silent(
#   y12_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                           target = ~ Y_12,
#                           data = sample_B2,
#                           pop_totals = X_totals[1:5],
#                           method_selection = "cloglog")
# )
# expect_equal(y12_corr_one$output$mean, 5.708058, tolerance = 0.0001) ## true value for this sim
# expect_equal(y12_corr_one$output$SE, 0.2604017, tolerance = 0.0001) ## true value for this sim
# expect_false(y12_corr_one$confidence_interval$lower_bound < mean(Y_12) &
#                y12_corr_one$confidence_interval$upper_bound > mean(Y_12)) ## conf int
#
#
# ## for y21
# expect_silent(
#   y21_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                           target = ~ Y_21,
#                           data = sample_B2,
#                           pop_totals = X_totals[1:5],
#                           method_selection = "cloglog")
# )
#
# expect_equal(y21_corr_one$output$mean, 0.6096105, tolerance = 0.0001) ## true value for this sim
# expect_equal(y21_corr_one$output$SE, 0.02698877, tolerance = 0.0001) ## true value for this sim
# expect_true(y21_corr_one$confidence_interval$lower_bound < mean(Y_21) &
#               y21_corr_one$confidence_interval$upper_bound > mean(Y_21)) ## conf int
#
# ## for y22
# expect_silent(
#   y22_corr_one <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                           target = ~ Y_22,
#                           data = sample_B2,
#                           pop_totals = X_totals[1:5],
#                           method_selection = "cloglog")
# )
#
# expect_equal(y22_corr_one$output$mean, 0.6210885, tolerance = 0.0001) ## true value for this sim
# expect_equal(y22_corr_one$output$SE, 0.02561999, tolerance = 0.0001) ## true value for this sim
# expect_true(y22_corr_one$confidence_interval$lower_bound < mean(Y_22) &
#               y22_corr_one$confidence_interval$upper_bound > mean(Y_22)) ## conf int
#
# ##### all target variables  ---------------------------------------------------------------
#
# expect_silent(
#   y_all_corr <- nonprob(selection = ~ X1 + X2 + X3 + X4,
#                         target = ~ Y_11 + Y_12 + Y_21 + Y_22,
#                         data = sample_B2,
#                         pop_totals = X_totals[1:5],
#                         method_selection = "cloglog")
# )
#
# expect_identical(y_all_corr$output$mean,
#                  c(y11_corr_one$output$mean, y12_corr_one$output$mean,
#                    y21_corr_one$output$mean, y22_corr_one$output$mean))
#
# expect_identical(y_all_corr$output$SE,
#                  c(y11_corr_one$output$SE, y12_corr_one$output$SE,
#                    y21_corr_one$output$SE, y22_corr_one$output$SE))
#
# expect_identical(y_all_corr$confidence_interval,
#                  data.frame(lower_bound = c(y11_corr_one$confidence_interval$lower_bound,
#                                             y12_corr_one$confidence_interval$lower_bound,
#                                             y21_corr_one$confidence_interval$lower_bound,
#                                             y22_corr_one$confidence_interval$lower_bound),
#                             upper_bound = c(y11_corr_one$confidence_interval$upper_bound,
#                                             y12_corr_one$confidence_interval$upper_bound,
#                                             y21_corr_one$confidence_interval$upper_bound,
#                                             y22_corr_one$confidence_interval$upper_bound),
#                             row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))





#### all X variables variables ------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

# ## for y11
# expect_silent(
#   y11_corr_all <- nonprob(selection = X_formula,
#                           target = ~ Y_11,
#                           data = sample_B2,
#                           pop_totals = X_totals,
#                           method_selection = "cloglog")
# )
#
# expect_equal(y11_corr_all$output$mean, 2.045187, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_all$output$SE, 0.1403925, tolerance = 0.0001) ## true value for this sim
# expect_true(y11_corr_all$confidence_interval$lower_bound < mean(Y_11) &
#               y11_corr_all$confidence_interval$upper_bound > mean(Y_11)) ## conf int
#
#
# ## for y12
# expect_silent(
#   y12_corr_all <- nonprob(selection = X_formula,
#                           target = ~ Y_12,
#                           data = sample_B2,
#                           pop_totals = X_totals,
#                           method_selection = "cloglog")
# )
# expect_equal(y12_corr_all$output$mean, 5.570961, tolerance = 0.0001) ## true value for this sim
# expect_equal(y12_corr_all$output$SE, 0.267995, tolerance = 0.0001) ## true value for this sim
# expect_false(y12_corr_all$confidence_interval$lower_bound < mean(Y_12) &
#                y12_corr_all$confidence_interval$upper_bound > mean(Y_12)) ## conf int
#
#
# ## for y21
# expect_silent(
#   y21_corr_all <- nonprob(selection = X_formula,
#                           target = ~ Y_21,
#                           data = sample_B2,
#                           pop_totals = X_totals,
#                           method_selection = "cloglog")
# )
#
# expect_equal(y21_corr_all$output$mean, 0.5999782, tolerance = 0.0001) ## true value for this sim
# expect_equal(y21_corr_all$output$SE, 0.02837228, tolerance = 0.0001) ## true value for this sim
# expect_true(y21_corr_all$confidence_interval$lower_bound < mean(Y_21) &
#               y21_corr_all$confidence_interval$upper_bound > mean(Y_21)) ## conf int
#
# ## for y22
# expect_silent(
#   y22_corr_all <- nonprob(selection = X_formula,
#                           target = ~ Y_22,
#                           data = sample_B2,
#                           pop_totals = X_totals,
#                           method_selection = "cloglog")
# )
#
# expect_equal(y22_corr_all$output$mean, 0.558637, tolerance = 0.0001) ## true value for this sim
# expect_equal(y22_corr_all$output$SE, 0.02463264, tolerance = 0.0001) ## true value for this sim
# expect_false(y22_corr_all$confidence_interval$lower_bound < mean(Y_22) &
#                y22_corr_all$confidence_interval$upper_bound > mean(Y_22)) ## conf int
#
# ##### all target variables  ---------------------------------------------------------------
#
# expect_silent(
#   y_all_corr_all <- nonprob(selection = X_formula,
#                             target = ~ Y_11 + Y_12 + Y_21 + Y_22,
#                             data = sample_B2,
#                             pop_totals = X_totals,
#                             method_selection = "cloglog")
# )
#
# expect_identical(y_all_corr_all$output$mean,
#                  c(y11_corr_all$output$mean, y12_corr_all$output$mean, y21_corr_all$output$mean, y22_corr_all$output$mean))
# #
# expect_identical(y_all_corr_all$output$SE,
#                  c(y11_corr_all$output$SE, y12_corr_all$output$SE, y21_corr_all$output$SE, y22_corr_all$output$SE))
# #
# expect_identical(y_all_corr_all$confidence_interval,
#                  data.frame(lower_bound = c(y11_corr_all$confidence_interval$lower_bound,
#                                             y12_corr_all$confidence_interval$lower_bound,
#                                             y21_corr_all$confidence_interval$lower_bound,
#                                             y22_corr_all$confidence_interval$lower_bound),
#                             upper_bound = c(y11_corr_all$confidence_interval$upper_bound,
#                                             y12_corr_all$confidence_interval$upper_bound,
#                                             y21_corr_all$confidence_interval$upper_bound,
#                                             y22_corr_all$confidence_interval$upper_bound),
#                             row.names = c("Y_11", "Y_12", "Y_21", "Y_22")))



#### variable selection  ------------------------------------------------------------------
##### one target variable  ----------------------------------------------------------------

## y_11
# expect_silent(
#   y11_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_11,
#                            data = sample_B2,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "cloglog",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "SCAD", nfolds = 5))
# )
#
# expect_equal(y11_corr_scad$output$mean, 1.992688, tolerance = 0.0001) ## true value for this sim
# expect_equal(y11_corr_scad$output$SE, 0.06352848, tolerance = 0.0001) ## true value for this sim
# expect_true(y11_corr_scad$confidence_interval$lower_bound < mean(Y_11) &
#               y11_corr_scad$confidence_interval$upper_bound > mean(Y_11)) ## conf int
# expect_true(NROW(y11_corr_scad$selection$coefficients) == 2)
#
# ## y_12
# expect_silent(
#   y12_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_12,
#                            data = sample_B2,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "cloglog",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "SCAD", nfolds = 5))
# )
#
# expect_equal(y12_corr_scad$output$mean, 5.712705, tolerance = 0.0001) ## true value for this sim
# expect_equal(y12_corr_scad$output$SE, 0.1460298, tolerance = 0.0001) ## true value for this sim
# expect_false(y12_corr_scad$confidence_interval$lower_bound < mean(Y_12) &
#                y12_corr_scad$confidence_interval$upper_bound > mean(Y_12)) ## conf int
# expect_true(NROW(y12_corr_scad$selection$coefficients) == 2)
#
# ## y_21
# expect_silent(
#   y21_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_21,
#                            data = sample_B2,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "cloglog",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "SCAD", nfolds = 5))
# )
#
# expect_equal(y21_corr_scad$output$mean, 0.5955036, tolerance = 0.0001) ## true value for this sim
# expect_equal(y21_corr_scad$output$SE, 0.01039547, tolerance = 0.0001) ## true value for this sim
# expect_false(y21_corr_scad$confidence_interval$lower_bound < mean(Y_21) &
#                y21_corr_scad$confidence_interval$upper_bound > mean(Y_21)) ## conf int
# expect_true(NROW(y21_corr_scad$selection$coefficients) == 2)
#
# ## y_22
# expect_silent(
#   y22_corr_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#                            target = ~ Y_22,
#                            data = sample_B2,
#                            pop_totals = X_totals[1:11],
#                            method_selection = "cloglog",
#                            control_inference = controlInf(vars_selection = TRUE),
#                            control_selection = controlSel(penalty = "SCAD", nfolds = 5))
# )
#
# expect_equal(y22_corr_scad$output$mean, 0.621987, tolerance = 0.0001) ## true value for this sim
# expect_equal(y22_corr_scad$output$SE, 0.01026539, tolerance = 0.0001) ## true value for this sim
# # to fix
# # expect_true(y22_corr_scad$confidence_interval$lower_bound < mean(Y_22) &
# #               y22_corr_scad$confidence_interval$upper_bound > mean(Y_22)) ## conf int
# expect_true(NROW(y22_corr_scad$selection$coefficients) == 2)
