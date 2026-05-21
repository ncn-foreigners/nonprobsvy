# nonprobsvy — developer conventions for Claude

## Sample labelling convention (IMPORTANT)

This package **always** uses the notation of the vignette and JSS paper:

- **Non-probability sample `S_{NP}`** — observes the outcome `y` and
  covariates `X`. Internally: selection indicator `R == 1`, suffix
  `_nons`, estimated propensity/selection score `ps`,
  inverse-probability weights `ipw_weights = 1/ps`.
- **Probability (reference) sample `S_{P}`** — observes only `X` (its
  `y` is never used), carries design weights. Internally: `R == 0`,
  suffix `_rand`, the `svydesign` object, `weights_rand`.

This matches Chen, Li & Wu (2020). The README, man pages
(`R/nonprob_documentation.R`, `R/method_*.R`, `R/control_outcome.R`) and
the vignette all use `S_{NP}` / `S_{P}`. Earlier revisions used `S_A`
(non-probability) / `S_B` (probability), so read **A = NP** and **B =
P** in older git history and commit messages.

> ⚠️ **Watch out — two key papers use the OPPOSITE labelling** (their
> **Sample A = probability**, **Sample B = non-probability**), so
> **translate the roles** when porting their equations into the code: -
> **Yang, Kim & Song (2020, JRSS-B)** — basis for the doubly-robust
> `bias_correction` path (`I_B`, `π_B` refer to the *non-probability*
> sampling score; joint estimating equations (9), variance (24)–(25)). -
> **Kim et al. (2021, JRSS-A)** — basis for the GLM mass-imputation
> analytic variance. In their `ĉ` (eq. (9)/(17)) the bread sum
> `n_B^{-1} Σ_{i∈B} ṁ h'` is over their *non-probability* sample (= our
> `S_{NP}`) and the meat `N^{-1} Σ_{i∈A} w_i ṁ` is over their
> *probability* sample (= our `S_{P}`).
>
> A mismatch here looks like a “wrong term in the equation” but is
> really just notation.

## Weights: keep two distinct concepts separate

- `case_weights` are **frequency weights** (how many population units a
  non-probability record represents). They are **NOT** propensity /
  inverse-probability weights.
- Inverse-probability (IP) weights are `1/ps` (= `ipw_weights`, or
  `bias_corr_ipw_weights` on the joint path).
- Do not substitute one for the other in estimating equations or
  variance formulas. In particular, the Yang–Kim–Song variance (eq. 25)
  uses `1/π_B`, **not** `case_weights`.
- Combining frequency weights *together with* IP weights is planned
  future work and is tracked in the issue tracker — do not improvise it
  inside variance code.
