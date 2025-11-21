# residualDiagnostics methods for [`lme`](https://rdrr.io/pkg/nlme/man/lme.html) objects

residualDiagnostics methods for
[`lme`](https://rdrr.io/pkg/nlme/man/lme.html) objects

## Usage

``` r
# S3 method for class 'lme'
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

  An object from [`lme`](https://rdrr.io/pkg/nlme/man/lme.html).

- ev.perc:

  The extreme value percentile to use. Defaults to .001.

- robust:

  A logical value, whether to use robust estimates or not. Defaults to
  `FALSE`.

- distr:

  A character string specifying the assumed distribution. Currently
  “normal”, but future options may be supported in the future.

- standardized:

  A logical value whether to use standardized pearson residual values or
  not. Defaults to `TRUE`.

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
m <- nlme::lme(extra ~ group, data = sleep, random = ~ 1 | ID,
  na.action = na.omit)

 residualDiagnostics(m)$Residuals
#>        Residuals  Predicted   isEV Index
#>            <num>      <num> <fctr> <int>
#>  1: -0.700128686 -0.9621668     No     2
#>  2:  0.021691421 -0.2197614     No     3
#>  3: -0.117253346 -1.0931795     No     4
#>  4:  0.658755591 -0.7001414     No     5
#>  5:  0.665700971  2.7935312     No     6
#>  6:  0.323895359  3.4049239     No     7
#>  7:  0.400316839  0.4353022     No     8
#>  8: -1.532411926  1.3960622     No     9
#>  9:  0.279433777  1.7454295     No    10
#> 10: -0.106066416  1.9966289     No    11
#> 11:  0.165397761  0.6493189     No    12
#> 12: -0.320216242  1.3917243     No    13
#> 13: -0.459161009  0.5183061     No    14
#> 14: -1.110119656  0.9113443     No    15
#> 15: -0.005506903  4.4050169     No    16
#> 16:  0.530821382  5.0164096     No    17
#> 17: -0.490424510  2.0467879     No    18
#> 18:  1.747982740  3.0075479     No    19
#> 19:  0.047292851  3.3569151     No    20

m <- nlme::lme(hp ~ mpg, data = mtcars, random = ~ 1 | cyl,
  na.action = na.omit)
residualDiagnostics(m)$Residuals
#>       Residuals Predicted   isEV Index
#>           <num>     <num> <fctr> <int>
#>  1: -0.23664029 118.76861     No     1
#>  2: -0.23664029 118.76861     No     2
#>  3: -0.22511638 101.34160     No     3
#>  4: -0.18972936 117.03034     No     4
#>  5: -0.43428124 191.09211     No     5
#>  6: -0.71168060 131.37103     No     6
#>  7:  0.93880308 210.21303     No     7
#>  8: -0.87407610  94.38853     No     8
#>  9: -0.17114196 101.34160     No     9
#> 10: -0.09690578 126.59080     No    10
#> 11: -0.26109405 132.67473     No    11
#> 12: -0.56908306 201.08714     No    12
#> 13: -0.46353346 197.17604     No    13
#> 14: -0.70981586 206.30194     No    14
#> 15: -0.59806685 227.16112     No    15
#> 16: -0.32819477 227.16112     No    16
#> 17:  0.58090589 208.47477     No    17
#> 18:  0.17209140  59.62322     No    18
#> 19: -0.44028418  68.31455     No    19
#> 20:  0.32102019  53.10473     No    20
#> 21: -0.26962808 106.99096     No    21
#> 22: -1.48424891 204.99824     No    22
#> 23: -1.51943211 206.30194     No    23
#> 24:  0.82152575 214.55870     No    24
#> 25: -0.37564257 188.91928     No    25
#> 26: -0.42602300  81.78611     No    26
#> 27:  0.09619667  87.43547     No    27
#> 28:  1.20593552  68.31455     No    28
#> 29:  1.62747604 203.69454     No    29
#> 30:  1.36506772 124.41797     No    30
#> 31:  3.44974596 207.17107    Yes    31
#> 32:  0.04249069 107.42552     No    32
#>       Residuals Predicted   isEV Index

rm(m, sleep)
```
