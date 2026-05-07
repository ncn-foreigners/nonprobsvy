# Admin Data (Non-Probability Survey)

This is a subset of the Central Job Offers Database, a voluntary
administrative data set (non-probability sample). The data was slightly
manipulated to ensure the relationships were preserved, and then
aligned. For more information about the CBOP, please refer to:
<https://oferty.praca.gov.pl/portal>.

## Usage

``` r
admin
```

## Format

A single data.frame with 9,344 rows and 6 columns

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

- `single_shift`:

  Whether an entity seeks employees on a single shift.

## Examples

``` r

data("admin")
head(admin)
#>    id private size nace region single_shift
#> 1 j_1       0    L    P     30        FALSE
#> 2 j_2       0    L    O     14         TRUE
#> 3 j_3       0    L    O     04         TRUE
#> 4 j_4       0    L    O     24         TRUE
#> 5 j_5       0    L    O     04         TRUE
#> 6 j_6       1    L    C     28        FALSE
```
