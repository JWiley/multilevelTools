# Calculate multilevel omega reliability

This function uses multilevel structural equation modelling to calculate
between and within reliability using coefficient omega.

## Usage

``` r
omegaSEM(items, id, data, savemodel = FALSE)
```

## Arguments

- items:

  A character vector giving the variables that map to the items in the
  scale. Note that these should be reverse scored prior to running this
  function.

- id:

  A character string giving the name of the variable that indicates
  which rows of the dataset belong to the same person or group for the
  multilevel analysis.

- data:

  A `data.table` or `data.frame` to be used for analysis.

- savemodel:

  A logical value indicating whether the underlying model should be
  saved and returned. Defaults to `FALSE`.

## Value

a list with two elements, the first, “Results” contains the estimates
for coefficient omega at the within and between level. The next element,
“Fit” contains the entire fitted model from lavaan, if
`savemodel = TRUE`.

## References

Geldhof, G. J., Preacher, K. J., & Zyphur, M. J. (2014)
\<doi:10.1037/a0032138\> "Reliability estimation in a multilevel
confirmatory factor analysis framework"

## Examples

``` r
  data(aces_daily, package = "JWileymisc")
  omegaSEM(
    items = c("COPEPrb", "COPEPrc", "COPEExp"),
    id = "UserID",
    data = aces_daily,
    savemodel = FALSE)
#> $Results
#>            label   est ci.lower ci.upper
#> 25  omega_within 0.719    0.697    0.740
#> 28 omega_between 0.908    0.884    0.933
#> 
```
