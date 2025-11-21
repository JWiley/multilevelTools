# residualDiagnostics methods for merMod objects

residualDiagnostics methods for merMod objects

## Usage

``` r
# S3 method for class 'merMod'
residualDiagnostics(
  object,
  ev.perc = 0.001,
  robust = FALSE,
  distr = "normal",
  standardized = TRUE,
  cut = 8L,
  quantiles = TRUE,
  ...
)
```

## Arguments

- object:

  An object with class
  [`merMod-class`](https://rdrr.io/pkg/lme4/man/merMod-class.html).
  Currently only [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) models
  are supported.

- ev.perc:

  The extreme value percentile to use. Defaults to .001.

- robust:

  A logical value, whether to use robust estimates or not. Defaults to
  `FALSE`.

- distr:

  A character string specifying the assumed distribution. Currently
  “normal”, but may expand in the future if
  [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) models are
  supported.

- standardized:

  A logical value whether to use standardized residual values or not.
  Defaults to `TRUE`.

- cut:

  An integer, how many unique predicted values there have to be at least
  for predicted values to be treated continuously, otherwise they are
  treated as discrete values. Defaults to 8.

- quantiles:

  A logical whether to calculate quantiles for the residuals. Defaults
  to `TRUE`. If `FALSE`, then do not calculate them. These are based on
  simple quantiles for each predicted value if the predicted values are
  few enough to be treated discretely. See `cut` argument. Otherwise
  they are based on quantile regression. First trying smoothing splines,
  and falling back to linear quantil regression if the splines fail. You
  may also want to turn these off if they are not working well, or are
  not of value in your diagnostics.

- ...:

  Additional arguments. Not currently used.

## Value

A logical
([`is.residualDiagnostics`](https://joshuawiley.com/JWileymisc/reference/residualDiagnostics.html))
or a residualDiagnostics object (list) for
[`as.residualDiagnostics`](https://joshuawiley.com/JWileymisc/reference/residualDiagnostics.html)
and
[`residualDiagnostics`](https://joshuawiley.com/JWileymisc/reference/residualDiagnostics.html).

## Examples

``` r
library(JWileymisc)
sleep[1,1] <- NA
m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep)

residualDiagnostics(m)$Residuals
#>        Residuals  Predicted   isEV Index
#>            <num>      <num> <fctr> <int>
#>  1: -0.700128956 -0.9621661     No     2
#>  2:  0.021690939 -0.2197610     No     3
#>  3: -0.117254135 -1.0931787     No     4
#>  4:  0.658754382 -0.7001407     No     5
#>  5:  0.665701421  2.7935304     No     6
#>  6:  0.323896367  3.4049228     No     7
#>  7:  0.400316373  0.4353024     No     8
#>  8: -1.532410421  1.3960619     No     9
#>  9:  0.279434030  1.7454290     No    10
#> 10: -0.106066615  1.9966292     No    11
#> 11:  0.165396845  0.6493196     No    12
#> 12: -0.320216425  1.3917247     No    13
#> 13: -0.459161498  0.5183069     No    14
#> 14: -1.110119449  0.9113449     No    15
#> 15: -0.005505897  4.4050160     No    16
#> 16:  0.530822260  5.0164084     No    17
#> 17: -0.490424247  2.0467880     No    18
#> 18:  1.747981709  3.0075476     No    19
#> 19:  0.047293317  3.3569147     No    20

#  gm1 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
#    data = lme4::cbpp, family = binomial)
# residualDiagnostics(gm1) ## should be an error

rm(m, sleep)
```
