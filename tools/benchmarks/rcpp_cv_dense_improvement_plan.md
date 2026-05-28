# Dense CV Improvement Plan

This note records the next optimization ideas for the `cpp-improvement` branch.
The branch currently contains only the dense C++ CV-fold parallelization in
`src/nonprobCV_cpp.cpp`.

## Current branch state

- Dense `arma::mat` implementation is preserved.
- CV folds are evaluated in parallel when `verbose = FALSE`.
- The current change preserves same-seed output versus the restored dense
  baseline: selected variables, selected lambda, final theta, and CV error
  vector.
- Earlier sparse/dense abstraction work was removed from this branch and should
  be treated as separate future work.

## Warm starts

`glmnet` fits a decreasing lambda path with warm starts: the previous lambda fit
is reused as the starting point for the next lambda. That idea was tested here
inside each CV fold.

Result: warm starts were not kept. They preserved selected variables, selected
lambda, and final theta in the trial cases, but changed the reported CV-error
vector. Under the strict "same results" rule, this fails the acceptance gate.

## Next experiments

1. Path-level allocation reuse

   Implement a fold-local path routine that still cold-starts every lambda fit,
   but reuses allocated matrices and vectors across lambda values. This should
   reduce allocation overhead without changing numerical starting points.

   Acceptance gate:
   - exact same selected variables, lambda, theta, and CV-error vector versus
     current `cpp-improvement`;
   - measurable speedup on a 50-variable dense simulation;
   - no R API change.

2. KKT-checked screening or active-set fitting

   Try a glmnet-inspired candidate-set strategy: fit on a reduced predictor set,
   compute the full score/check over all predictors, and refit with the full set
   if an excluded variable violates the check.

   Acceptance gate:
   - exact same selected variables, lambda, theta, and CV-error vector;
   - fallback must guarantee equivalence when screening is unsafe;
   - likely useful only when the number of predictors is larger than 50.

3. Avoid repeated fold matrix copies

   Investigate whether training/test row views or precomputed fold matrices can
   reduce copying without hurting dense BLAS performance.

   Acceptance gate:
   - exact same results;
   - benchmark must show a real gain, not just less copying in theory.

## Ideas to avoid for now

- Warm starts: rejected for this branch because they changed the CV-error vector.
- Early stopping along the lambda path: likely changes reported diagnostics.
- Full coordinate-descent rewrite: too invasive and unlikely to preserve current
  SCAD/MCP/lasso estimating-equation behavior exactly.

