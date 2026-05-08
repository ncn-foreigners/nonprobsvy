# Rcpp CV Benchmark Report

- Timestamp: `20260508-213853`
- Baseline ref: `origin/dev`
- First-pass ref: `194373e`
- Current ref: worktree at `/Users/berenz/mac/nauka/ncn-foreigners/software/nonprobsvy`
- Repetitions per scenario: `5`
- Quick mode: `FALSE`
- Required median direct-CV speedup versus dev: `2.00x`
- Incremental speed threshold versus first pass: `1.10x`
- Incremental max-RSS reduction threshold versus first pass: `10.0%`

## Gate

- Median direct-CV speedup versus dev: `7.575x`
- Median direct-CV speedup versus first pass: `1.435x`
- Median direct-CV max-RSS reduction versus first pass: `-0.2%`
- Numerical equivalence versus first pass (`max_abs_diff <= 1e-6` and identical selected sets): `TRUE`
- Numerical equivalence versus dev (`max_abs_diff <= 1e-6` and identical selected sets): `TRUE`
- End-to-end regression-free versus first pass (`speedup >= 0.95x`): `TRUE`
- Incremental runtime pass: `TRUE`
- Incremental memory pass: `FALSE`
- Overall incremental gate: `PASS`

## Scenario Results

- `direct_cv|logit|moderate|pop=FALSE|est=none`: time dev `0.1220s`, first `0.0280s`, current `0.0170s`; speedup current/dev `7.176x`, current/first `1.647x`; RSS dev `276.6 MB`, first `240.4 MB`, current `239.2 MB`; RSS reduction current/first `0.5%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `7.841e-16`, selected vs first `TRUE`
- `direct_cv|logit|moderate|pop=TRUE|est=none`: time dev `0.1100s`, first `0.0240s`, current `0.0150s`; speedup current/dev `7.333x`, current/first `1.600x`; RSS dev `274.9 MB`, first `233.9 MB`, current `233.8 MB`; RSS reduction current/first `0.0%`; max fitted diff vs first `5.204e-18`, CV-error diff vs first `2.571e-08`, selected vs first `TRUE`
- `direct_cv|logit|high|pop=FALSE|est=none`: time dev `0.3180s`, first `0.0460s`, current `0.0350s`; speedup current/dev `9.086x`, current/first `1.314x`; RSS dev `271.0 MB`, first `233.2 MB`, current `234.2 MB`; RSS reduction current/first `-0.4%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `2.961e-03`, selected vs first `TRUE`
- `direct_cv|logit|high|pop=TRUE|est=none`: time dev `0.4010s`, first `0.0590s`, current `0.0430s`; speedup current/dev `9.326x`, current/first `1.372x`; RSS dev `265.3 MB`, first `234.4 MB`, current `235.6 MB`; RSS reduction current/first `-0.5%`; max fitted diff vs first `1.110e-16`, CV-error diff vs first `1.055e-15`, selected vs first `TRUE`
- `end_to_end|logit|moderate|pop=FALSE|est=mle`: time dev `0.0450s`, first `0.0280s`, current `0.0280s`; speedup current/dev `1.607x`, current/first `1.000x`; RSS dev `254.0 MB`, first `247.3 MB`, current `246.6 MB`; RSS reduction current/first `0.3%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `NA`, selected vs first `TRUE`
- `end_to_end|logit|moderate|pop=FALSE|est=gee`: time dev `0.0280s`, first `0.0110s`, current `0.0100s`; speedup current/dev `2.800x`, current/first `1.100x`; RSS dev `244.1 MB`, first `238.2 MB`, current `242.5 MB`; RSS reduction current/first `-1.8%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `NA`, selected vs first `TRUE`
- `direct_cv|probit|moderate|pop=FALSE|est=none`: time dev `0.0800s`, first `0.0220s`, current `0.0140s`; speedup current/dev `5.714x`, current/first `1.571x`; RSS dev `260.0 MB`, first `235.1 MB`, current `234.4 MB`; RSS reduction current/first `0.3%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `7.649e-05`, selected vs first `TRUE`
- `direct_cv|probit|moderate|pop=TRUE|est=none`: time dev `0.0980s`, first `0.0250s`, current `0.0170s`; speedup current/dev `5.765x`, current/first `1.471x`; RSS dev `258.1 MB`, first `235.0 MB`, current `231.3 MB`; RSS reduction current/first `1.6%`; max fitted diff vs first `3.469e-18`, CV-error diff vs first `1.826e-12`, selected vs first `TRUE`
- `direct_cv|probit|high|pop=FALSE|est=none`: time dev `0.2080s`, first `0.0380s`, current `0.0280s`; speedup current/dev `7.429x`, current/first `1.357x`; RSS dev `261.3 MB`, first `234.4 MB`, current `235.3 MB`; RSS reduction current/first `-0.4%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `2.552e-07`, selected vs first `TRUE`
- `direct_cv|probit|high|pop=TRUE|est=none`: time dev `0.3580s`, first `0.0590s`, current `0.0460s`; speedup current/dev `7.783x`, current/first `1.283x`; RSS dev `294.5 MB`, first `233.4 MB`, current `237.4 MB`; RSS reduction current/first `-1.7%`; max fitted diff vs first `4.163e-17`, CV-error diff vs first `4.957e+00`, selected vs first `TRUE`
- `end_to_end|probit|moderate|pop=FALSE|est=mle`: time dev `0.0430s`, first `0.0320s`, current `0.0310s`; speedup current/dev `1.387x`, current/first `1.032x`; RSS dev `270.4 MB`, first `254.9 MB`, current `260.6 MB`; RSS reduction current/first `-2.2%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `NA`, selected vs first `TRUE`
- `end_to_end|probit|moderate|pop=FALSE|est=gee`: time dev `0.0250s`, first `0.0130s`, current `0.0120s`; speedup current/dev `2.083x`, current/first `1.083x`; RSS dev `245.2 MB`, first `241.5 MB`, current `244.7 MB`; RSS reduction current/first `-1.3%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `NA`, selected vs first `TRUE`
- `direct_cv|cloglog|moderate|pop=FALSE|est=none`: time dev `0.1550s`, first `0.0310s`, current `0.0210s`; speedup current/dev `7.381x`, current/first `1.476x`; RSS dev `283.6 MB`, first `232.8 MB`, current `238.2 MB`; RSS reduction current/first `-2.3%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `8.237e-04`, selected vs first `TRUE`
- `direct_cv|cloglog|moderate|pop=TRUE|est=none`: time dev `0.1390s`, first `0.0270s`, current `0.0180s`; speedup current/dev `7.722x`, current/first `1.500x`; RSS dev `272.4 MB`, first `233.5 MB`, current `237.5 MB`; RSS reduction current/first `-1.7%`; max fitted diff vs first `1.735e-18`, CV-error diff vs first `2.940e-13`, selected vs first `TRUE`
- `direct_cv|cloglog|high|pop=FALSE|est=none`: time dev `0.3600s`, first `0.0580s`, current `0.0430s`; speedup current/dev `8.372x`, current/first `1.349x`; RSS dev `297.8 MB`, first `234.4 MB`, current `232.6 MB`; RSS reduction current/first `0.8%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `5.139e-03`, selected vs first `TRUE`
- `direct_cv|cloglog|high|pop=TRUE|est=none`: time dev `0.4400s`, first `0.0700s`, current `0.0500s`; speedup current/dev `8.800x`, current/first `1.400x`; RSS dev `311.8 MB`, first `238.1 MB`, current `238.0 MB`; RSS reduction current/first `0.0%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `4.035e-08`, selected vs first `TRUE`
- `end_to_end|cloglog|moderate|pop=FALSE|est=mle`: time dev `0.0340s`, first `0.0110s`, current `0.0110s`; speedup current/dev `3.091x`, current/first `1.000x`; RSS dev `251.0 MB`, first `240.9 MB`, current `235.2 MB`; RSS reduction current/first `2.4%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `NA`, selected vs first `TRUE`
- `end_to_end|cloglog|moderate|pop=FALSE|est=gee`: time dev `0.0360s`, first `0.0130s`, current `0.0110s`; speedup current/dev `3.273x`, current/first `1.182x`; RSS dev `251.1 MB`, first `236.0 MB`, current `240.0 MB`; RSS reduction current/first `-1.7%`; max fitted diff vs first `0.000e+00`, CV-error diff vs first `NA`, selected vs first `TRUE`

## Artifacts

- CSV: `rcpp_cv_benchmark_20260508-213853.csv`
- Latest CSV copy: `rcpp_cv_benchmark_latest.csv`
