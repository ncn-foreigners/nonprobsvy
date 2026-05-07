source("_code_for_all_.R")

# standard DR analytic variance should agree with one-outcome fits ----------

dr_glm_y1 <- nonprob(
  selection = ~x1 + x2,
  outcome = y1 ~ x1 + x2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  pop_size = sum(weights(kim2019_sample_prob)),
  method_selection = "logit",
  method_outcome = "glm"
)

dr_glm_y2 <- nonprob(
  selection = ~x1 + x2,
  outcome = y2 ~ x1 + x2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  pop_size = sum(weights(kim2019_sample_prob)),
  method_selection = "logit",
  method_outcome = "glm"
)

dr_glm_joint <- nonprob(
  selection = ~x1 + x2,
  outcome = y1 + y2 ~ x1 + x2,
  svydesign = kim2019_sample_prob,
  data = kim2019_sample_nonprob,
  pop_size = sum(weights(kim2019_sample_prob)),
  method_selection = "logit",
  method_outcome = "glm"
)

expect_equal(
  dr_glm_joint$output["y1", , drop = FALSE],
  dr_glm_y1$output
)

expect_equal(
  dr_glm_joint$output["y2", , drop = FALSE],
  dr_glm_y2$output
)

expect_equal(
  dr_glm_joint$confidence_interval["y1", , drop = FALSE],
  dr_glm_y1$confidence_interval
)

expect_equal(
  dr_glm_joint$confidence_interval["y2", , drop = FALSE],
  dr_glm_y2$confidence_interval
)

# combined DR analytic variance should agree with one-outcome fits ----------

admin_multi <- transform(
  admin,
  y1 = as.numeric(single_shift),
  y2 = as.numeric(single_shift) + 0.25 * as.numeric(private)
)

dr_comb_inf <- control_inf(vars_selection = TRUE, vars_combine = TRUE)
dr_comb_out <- control_out(nfolds = 2, nlambda = 5)
dr_comb_sel <- control_sel(nfolds = 2, nlambda = 5)

set.seed(2024)
dr_comb_y1 <- nonprob(
  selection = ~nace + size,
  outcome = y1 ~ nace + size,
  svydesign = jvs_svy,
  data = admin_multi,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  method_outcome = "glm",
  family_outcome = "gaussian",
  control_inference = dr_comb_inf,
  control_outcome = dr_comb_out,
  control_selection = dr_comb_sel
)

set.seed(2024)
dr_comb_y2 <- nonprob(
  selection = ~nace + size,
  outcome = y2 ~ nace + size,
  svydesign = jvs_svy,
  data = admin_multi,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  method_outcome = "glm",
  family_outcome = "gaussian",
  control_inference = dr_comb_inf,
  control_outcome = dr_comb_out,
  control_selection = dr_comb_sel
)

set.seed(2024)
dr_comb_joint <- nonprob(
  selection = ~nace + size,
  outcome = y1 + y2 ~ nace + size,
  svydesign = jvs_svy,
  data = admin_multi,
  pop_size = sum(weights(jvs_svy)),
  method_selection = "logit",
  method_outcome = "glm",
  family_outcome = "gaussian",
  control_inference = dr_comb_inf,
  control_outcome = dr_comb_out,
  control_selection = dr_comb_sel
)

expect_equal(
  dr_comb_joint$output["y1", , drop = FALSE],
  dr_comb_y1$output,
  tolerance = 1e-8
)

expect_equal(
  dr_comb_joint$output["y2", , drop = FALSE],
  dr_comb_y2$output,
  tolerance = 1e-8
)

expect_equal(
  dr_comb_joint$confidence_interval["y1", , drop = FALSE],
  dr_comb_y1$confidence_interval,
  tolerance = 1e-8
)

expect_equal(
  dr_comb_joint$confidence_interval["y2", , drop = FALSE],
  dr_comb_y2$confidence_interval,
  tolerance = 1e-8
)
