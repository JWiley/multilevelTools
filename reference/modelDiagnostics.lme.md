# modelDiagnostics method for lme objects

This function creates a number of diagnostics for
[`lme`](https://rdrr.io/pkg/nlme/man/lme.html) models.

## Usage

``` r
# S3 method for class 'lme'
modelDiagnostics(
  object,
  ev.perc = 0.001,
  robust = FALSE,
  distr = "normal",
  standardized = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted model object from
  [`lme`](https://rdrr.io/pkg/nlme/man/lme.html).

- ev.perc:

  A real number between 0 and 1 indicating the proportion of the
  theoretical distribution beyond which values are considered extreme
  values (possible outliers). Defaults to .001.

- robust:

  Whether to use robust mean and standard deviation estimates for normal
  distribution

- distr:

  A character string given the assumed distribution. Passed on to
  [`testDistribution`](https://joshuawiley.com/JWileymisc/reference/testDistribution.html).
  Defaults to “normal”.

- standardized:

  A logical whether to use standardized pearson residuals. Defaults to
  `TRUE` generally where possible but may depend on method.

- ...:

  Additional arguments, passed to
  [`residualDiagnostics`](https://joshuawiley.com/JWileymisc/reference/residualDiagnostics.html).

## Value

A logical
([`is.modelDiagnostics`](https://joshuawiley.com/JWileymisc/reference/modelDiagnostics.html))
or a modelDiagnostics object (list) for
[`as.modelDiagnostics`](https://joshuawiley.com/JWileymisc/reference/modelDiagnostics.html)
and
[`modelDiagnostics`](https://joshuawiley.com/JWileymisc/reference/modelDiagnostics.html).

## Examples

``` r
library(JWileymisc)
sleep[1,1] <- NA
m <- nlme::lme(extra ~ group, data = sleep,
 random = ~ 1 | ID, na.action = "na.omit")

 md <- modelDiagnostics(m, ev.perc = .1)
 md$extremeValues
#>    extra     ID Index                     EffectType
#>    <num> <fctr> <int>                         <char>
#> 1:   0.0      9     9                      Residuals
#> 2:  -0.1      5    15                      Residuals
#> 3:   4.6      9    19                      Residuals
#> 4:   3.4      6     6 Random Effect ID : (Intercept)
#> 5:   3.7      7     7 Random Effect ID : (Intercept)
#> 6:   4.4      6    16 Random Effect ID : (Intercept)
#> 7:   5.5      7    17 Random Effect ID : (Intercept)

 plot(md)




data(aces_daily, package = "JWileymisc")
m <- nlme::lme(PosAff ~ STRESS, data = aces_daily,
  random = ~ 1 + STRESS | UserID, na.action = "na.omit")
md <- modelDiagnostics(m, ev.perc = .001)
md$extremeValues
#>       PosAff UserID Index                        EffectType
#>        <num>  <int> <int>                            <char>
#>  1: 2.018930     19   643                         Residuals
#>  2: 4.647537     22   736                         Residuals
#>  3: 4.379594     22   737                         Residuals
#>  4: 4.285682     51  1759                         Residuals
#>  5: 5.000000     52  1819                         Residuals
#>  6: 4.593965     53  1838                         Residuals
#>  7: 5.000000     53  1841                         Residuals
#>  8: 3.679191     56  1934                         Residuals
#>  9: 3.499317     69  2371                         Residuals
#> 10: 1.333848     75  2580                         Residuals
#> 11: 4.200560     78  2688                         Residuals
#> 12: 1.000000     81  2805                         Residuals
#> 13: 1.723817     83  2884                         Residuals
#> 14: 4.702814     88  3018                         Residuals
#> 15: 4.982255     97  3328                         Residuals
#> 16: 4.851732     97  3329                         Residuals
#> 17: 4.697991     97  3330                         Residuals
#> 18: 5.000000    107  3675                         Residuals
#> 19: 1.173068    113  3902                         Residuals
#> 20: 1.000000    115  3954                         Residuals
#> 21: 1.041260    119  4100                         Residuals
#> 22: 4.713139    127  4357                         Residuals
#> 23: 4.054051    141  4834                         Residuals
#> 24: 3.975565    141  4835                         Residuals
#> 25: 3.989508    141  4836                         Residuals
#> 26: 4.278470    141  4844                         Residuals
#> 27: 4.447151    147  5079                         Residuals
#> 28: 1.079654    155  5341                         Residuals
#> 29: 5.000000    155  5358                         Residuals
#> 30: 4.660743    160  5520                         Residuals
#> 31: 1.342503    163  5626                         Residuals
#> 32: 4.706810    173  5944                         Residuals
#> 33: 4.194950    180  6212                         Residuals
#> 34: 4.258350    191  6568                         Residuals
#> 35: 5.000000    191  6569                         Residuals
#> 36: 4.701488    123  4219 Multivariate Random Effect UserID
#> 37: 3.942647    123  4220 Multivariate Random Effect UserID
#> 38: 3.841534    123  4221 Multivariate Random Effect UserID
#> 39: 3.072598    123  4222 Multivariate Random Effect UserID
#> 40: 4.165724    123  4223 Multivariate Random Effect UserID
#> 41: 2.573955    123  4224 Multivariate Random Effect UserID
#> 42: 1.871370    123  4225 Multivariate Random Effect UserID
#> 43: 4.497994    123  4226 Multivariate Random Effect UserID
#> 44: 2.207134    123  4227 Multivariate Random Effect UserID
#> 45: 4.370643    123  4228 Multivariate Random Effect UserID
#> 46: 1.713198    123  4229 Multivariate Random Effect UserID
#> 47: 3.786528    123  4230 Multivariate Random Effect UserID
#> 48: 4.702243    123  4231 Multivariate Random Effect UserID
#> 49: 2.262491    123  4232 Multivariate Random Effect UserID
#> 50: 5.000000    123  4234 Multivariate Random Effect UserID
#> 51: 2.392975    123  4235 Multivariate Random Effect UserID
#> 52: 3.588502    123  4236 Multivariate Random Effect UserID
#> 53: 3.446311    123  4237 Multivariate Random Effect UserID
#> 54: 2.832955    123  4238 Multivariate Random Effect UserID
#> 55: 3.671293    123  4239 Multivariate Random Effect UserID
#> 56: 4.893790    123  4240 Multivariate Random Effect UserID
#> 57: 4.633870    123  4242 Multivariate Random Effect UserID
#>       PosAff UserID Index                        EffectType
plot(md)







m <- nlme::lme(extra ~ 1, data = sleep, random = ~ 1 | ID/group,
  na.action = "na.omit")

md <- modelDiagnostics(m, ev.perc = .1)
md$extremeValues
#>    extra     ID  group Index                     EffectType
#>    <num> <fctr> <fctr> <int>                         <char>
#> 1:  -1.6      2      1     2                      Residuals
#> 2:   0.0      9      1     9                      Residuals
#> 3:   5.5      7      2    17                      Residuals
#> 4:   4.6      9      2    19                      Residuals
#> 5:   3.4      6      1     6 Random Effect ID : (Intercept)
#> 6:   3.7      7      1     7 Random Effect ID : (Intercept)
#> 7:   4.4      6      2    16 Random Effect ID : (Intercept)
#> 8:   5.5      7      2    17 Random Effect ID : (Intercept)
plot(md)




rm(m, md, sleep)
```
