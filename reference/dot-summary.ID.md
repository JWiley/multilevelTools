# Summarize a Variable in a Long Dataset by ID

Summarize a Variable in a Long Dataset by ID

## Usage

``` r
.summary.ID(data, var, idvar, CI = 0.95, robust = FALSE)
```

## Arguments

- data:

  A data.table object, long format

- var:

  A character string, the name of the variable to summarize

- idvar:

  A character string, the name of the grouping variable

- CI:

  A numeric value, the confidence interval to use. Default is .95.

- robust:

  A logical. Default is `FALSE`. If `TRUE`, the function will use the
  median as the estimate. If `FALSE`, the function will use the mean as
  the estimate.

## Value

A data.table object with the mean/median, lower limit, and upper limit
of the variable specified in `var` for each level of the grouping
variable specified in `idvar`.
