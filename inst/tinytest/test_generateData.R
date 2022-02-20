## checking functions from generateData.r
## general
#### size
expect_equal(
  nrow(genSimData(N=10000,n=1000)), 10000
)

#### rho equal to n
expect_equal(
  sum(genSimData(N=10000,n=1000)$rho),
  1000,
  tolerance = 0.01
)

### relation max(srs)/min(srs) = 50

expect_equal(
  {
    z <- genSimData(N=10000,n=1000)
    max(z$srs)/min(z$srs)
  },
  50,
  tolerance = 0.01
)

## the same dataset with the same seed

expect_equal(
  {
    set.seed(123)
    genSimData(N=10000,n=1000)
  },
  {
    set.seed(123)
    genSimData(N=10000,n=1000)
  }
)
