# Intraclass Correlation Coefficient (ICC) from Mixed Models

This function estimates the ICC from mixed effects models estimated
using lme4.

## Usage

``` r
iccMixed(dv, id, data, family = c("gaussian", "binomial"))
```

## Arguments

- dv:

  A character string giving the variable name of the dependent variable.

- id:

  A character vector of length one or more giving the ID variable(s).
  Can be more than one.

- data:

  A data.table containing the variables used in the formula. This is a
  required argument. If a data.frame, it will silently coerce to a
  data.table. If not a data.table or data.frame, it will attempt to
  coerce, with a message.

- family:

  A character vector giving the family to use for the model. Currently
  only supports “gaussian” or “binomial”.

## Value

A data table of the ICCs

## References

For details, see Campbell, M. K., Mollison, J., and Grimshaw, J. M.
(2001) \<doi:10.1002/1097-0258(20010215)20:3 "Cluster trials in
implementation research: estimation of intracluster correlation
coefficients and sample size."

## Examples

``` r
iccMixed("mpg", "cyl", mtcars)
#>         Var    Sigma       ICC
#>      <char>    <num>     <num>
#> 1:      cyl 33.19401 0.7617095
#> 2: Residual 10.38430 0.2382905
iccMixed("mpg", "cyl", data.table::as.data.table(mtcars))
#>         Var    Sigma       ICC
#>      <char>    <num>     <num>
#> 1:      cyl 33.19401 0.7617095
#> 2: Residual 10.38430 0.2382905
iccMixed("mpg", "cyl", data.table::as.data.table(mtcars), family = "gaussian")
#>         Var    Sigma       ICC
#>      <char>    <num>     <num>
#> 1:      cyl 33.19401 0.7617095
#> 2: Residual 10.38430 0.2382905
iccMixed("mpg", c("cyl", "am"), data.table::as.data.table(mtcars))
#>         Var     Sigma        ICC
#>      <char>     <num>      <num>
#> 1:      cyl 26.532017 0.68284241
#> 2:       am  2.921697 0.07519438
#> 3: Residual  9.401542 0.24196321
iccMixed("am", "cyl", data.table::as.data.table(mtcars), family = "binomial")
#>         Var     Sigma       ICC
#>      <char>     <num>     <num>
#> 1:      cyl 0.9467022 0.2234596
#> 2: Residual 3.2898681 0.7765404
```
