# Job Vacancy Survey

This is a subset of the Job Vacancy Survey from Poland (for one
quarter). The data has been subject to slight manipulation, but the
relationships in the data have been preserved. For further details on
the JVS, please refer to the following link:
<https://stat.gov.pl/obszary-tematyczne/rynek-pracy/popyt-na-prace/zeszyt-metodologiczny-popyt-na-prace,3,1.html>.

## Usage

``` r
jvs
```

## Format

A single data.frame with 6,523 rows and 6 columns

- `id`:

  Identifier of an entity (company: legal or local).

- `private`:

  Whether the company is a private (1) or public (0) entity.

- `size`:

  The size of the entity: S – small (to 9 employees), M – medium (10-49)
  or L – large (over 49).

- `nace`:

  The main NACE code for a given entity: C, D.E, F, G, H, I, J, K.L, M,
  N, O, P, Q or R.S (14 levels, 3 combined: D and E, K and L, and R and
  S).

- `region`:

  The region of Poland (16 levels: 02, 04, ..., 32).

- `weight`:

  The final (calibrated) weight (w-weight). We do not have access to
  design weights (d-weights).

## Examples

``` r

data("jvs")
head(jvs)
#>    id private size nace region weight
#> 1 j_1       0    L    O     14      1
#> 2 j_2       0    L    O     24      6
#> 3 j_3       0    L  R.S     14      1
#> 4 j_4       0    L  R.S     14      1
#> 5 j_5       0    L  R.S     22      1
#> 6 j_6       0    M  R.S     26      1
```
