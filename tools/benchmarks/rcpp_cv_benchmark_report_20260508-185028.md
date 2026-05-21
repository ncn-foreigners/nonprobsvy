# Rcpp CV Benchmark Report

- Timestamp: `20260508-185028`
- Baseline ref: `origin/dev`
- Current ref: worktree at `/Users/berenz/mac/nauka/ncn-foreigners/software/nonprobsvy`
- Repetitions per scenario: `3`
- Quick mode: `FALSE`
- Required median direct-CV speedup: `2.00x`

## Gate

- Median direct-CV speedup: `5.191x`
- Numerical equivalence pass (`max_abs_diff <= 1e-6` for fitted outputs and identical selected sets): `TRUE`
- End-to-end regression-free pass (`speedup >= 0.95` for representative workloads): `TRUE`
- Overall gate: `PASS`

## Scenario Results

- `direct_cv|logit|moderate|pop=FALSE|est=none`: baseline `0.1220s`, current `0.0270s`, speedup `4.519x`, max fitted diff `0`, CV-error diff `5.551115e-16`, selected identical `TRUE`
- `direct_cv|logit|moderate|pop=TRUE|est=none`: baseline `0.1090s`, current `0.0240s`, speedup `4.542x`, max fitted diff `1.73472e-18`, CV-error diff `2.629334e-08`, selected identical `TRUE`
- `direct_cv|logit|high|pop=FALSE|est=none`: baseline `0.3170s`, current `0.0500s`, speedup `6.340x`, max fitted diff `0`, CV-error diff `2.768993e-04`, selected identical `TRUE`
- `direct_cv|logit|high|pop=TRUE|est=none`: baseline `0.3920s`, current `0.0650s`, speedup `6.031x`, max fitted diff `5.55112e-17`, CV-error diff `4.440892e-16`, selected identical `TRUE`
- `end_to_end|logit|moderate|pop=FALSE|est=mle`: baseline `0.0450s`, current `0.0280s`, speedup `1.607x`, max fitted diff `0`, CV-error diff `NA`, selected identical `TRUE`
- `end_to_end|logit|moderate|pop=FALSE|est=gee`: baseline `0.0280s`, current `0.0100s`, speedup `2.800x`, max fitted diff `0`, CV-error diff `NA`, selected identical `TRUE`
- `direct_cv|probit|moderate|pop=FALSE|est=none`: baseline `0.0800s`, current `0.0220s`, speedup `3.636x`, max fitted diff `0`, CV-error diff `5.906672e-05`, selected identical `TRUE`
- `direct_cv|probit|moderate|pop=TRUE|est=none`: baseline `0.0980s`, current `0.0250s`, speedup `3.920x`, max fitted diff `1.73472e-18`, CV-error diff `1.123490e-12`, selected identical `TRUE`
- `direct_cv|probit|high|pop=FALSE|est=none`: baseline `0.2070s`, current `0.0390s`, speedup `5.308x`, max fitted diff `0`, CV-error diff `1.790737e-07`, selected identical `TRUE`
- `direct_cv|probit|high|pop=TRUE|est=none`: baseline `0.3490s`, current `0.0600s`, speedup `5.817x`, max fitted diff `0`, CV-error diff `2.413010e-01`, selected identical `TRUE`
- `end_to_end|probit|moderate|pop=FALSE|est=mle`: baseline `0.0420s`, current `0.0320s`, speedup `1.312x`, max fitted diff `0`, CV-error diff `NA`, selected identical `TRUE`
- `end_to_end|probit|moderate|pop=FALSE|est=gee`: baseline `0.0240s`, current `0.0130s`, speedup `1.846x`, max fitted diff `0`, CV-error diff `NA`, selected identical `TRUE`
- `direct_cv|cloglog|moderate|pop=FALSE|est=none`: baseline `0.1560s`, current `0.0310s`, speedup `5.032x`, max fitted diff `0`, CV-error diff `1.528367e-03`, selected identical `TRUE`
- `direct_cv|cloglog|moderate|pop=TRUE|est=none`: baseline `0.1370s`, current `0.0270s`, speedup `5.074x`, max fitted diff `0`, CV-error diff `2.146616e-13`, selected identical `TRUE`
- `direct_cv|cloglog|high|pop=FALSE|est=none`: baseline `0.3570s`, current `0.0570s`, speedup `6.263x`, max fitted diff `0`, CV-error diff `5.857482e-04`, selected identical `TRUE`
- `direct_cv|cloglog|high|pop=TRUE|est=none`: baseline `0.4370s`, current `0.0740s`, speedup `5.905x`, max fitted diff `0`, CV-error diff `4.217423e-05`, selected identical `TRUE`
- `end_to_end|cloglog|moderate|pop=FALSE|est=mle`: baseline `0.0340s`, current `0.0110s`, speedup `3.091x`, max fitted diff `0`, CV-error diff `NA`, selected identical `TRUE`
- `end_to_end|cloglog|moderate|pop=FALSE|est=gee`: baseline `0.0350s`, current `0.0120s`, speedup `2.917x`, max fitted diff `0`, CV-error diff `NA`, selected identical `TRUE`

## Artifacts

- CSV: `rcpp_cv_benchmark_20260508-185028.csv`
- Latest CSV copy: `rcpp_cv_benchmark_latest.csv`
