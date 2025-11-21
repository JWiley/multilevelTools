# Convert ranef() output to long format

Convert ranef() output to long format

## Usage

``` r
ranef2long(x, idvar)

.re.data(d, i, idvar)
```

## Arguments

- x:

  A `brmsfit` object

- idvar:

  A character string specifying the grouping variable name for the
  random effects.

- d:

  A `ranef` object

- i:

  an integer, which random effect to pull out

## Value

A data.table object with the random effects in long format.
