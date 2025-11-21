# Estimate the effective sample size from longitudinal data

This function estimates the (approximate) effective sample size.

## Usage

``` r
nEffective(n, k, icc, dv, id, data, family = c("gaussian", "binomial"))
```

## Arguments

- n:

  The number of unique/independent units of observation

- k:

  The (average) number of observations per unit

- icc:

  The estimated ICC. If missing, will estimate (and requires that the
  family argument be correctly specified).

- dv:

  A character string giving the variable name of the dependent variable.

- id:

  A character vector of length one giving the ID variable.

- data:

  A data.table containing the variables used in the formula. This is a
  required argument. If a data.frame, it will silently coerce to a
  data.table. If not a data.table or data.frame, it will attempt to
  coerce, with a message.

- family:

  A character vector giving the family to use for the model. Currently
  only supports “gaussian” or “binomial”.

## Value

A data.table including the effective sample size.

## References

For details, see Campbell, M. K., Mollison, J., and Grimshaw, J. M.
(2001) \<doi:10.1002/1097-0258(20010215)20:3 "Cluster trials in
implementation research: estimation of intracluster correlation
coefficients and sample size."

## Examples

``` r
## example where n, k, and icc are estimated from the data
## provided, partly using iccMixed function
nEffective(dv = "mpg", id = "cyl", data = mtcars)
#>                     Type         N
#>                   <char>     <num>
#> 1: Effective Sample Size  3.826291
#> 2:     Independent Units  3.000000
#> 3:    Total Observations 32.000000

## example where n, k, and icc are known (or being 'set')
## useful for sensitivity analyses
nEffective(n = 60, k = 10, icc = .6)
#>                     Type      N
#>                   <char>  <num>
#> 1: Effective Sample Size  93.75
#> 2:     Independent Units  60.00
#> 3:    Total Observations 600.00
```
