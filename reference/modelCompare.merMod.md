# Compare two lmer models

This function provides fit statistics and effect sizes for model
comparisons. The models must be nested.

## Usage

``` r
# S3 method for class 'merMod'
modelCompare(model1, model2, ...)
```

## Arguments

- model1:

  A model estimated by [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

- model2:

  A model estimated by [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

- ...:

  Additional arguments, not currently used but included to match
  generic.

## Value

a data table with the fit indices for each model and comparing models to
each other.

## References

For estimating the marginal and conditional R-squared values, see:
Nakagawa, S. and Schielzeth, H. (2013). A general and simple method for
obtaining R2 from generalized linear mixed-effects models. Methods in
Ecology and Evolution, 4(2), 133-142. as well as: Johnson, P. C. (2014).
Extension of Nakagawa & Schielzeth's R2GLMM to random slopes models.
Methods in Ecology and Evolution, 5(9), 944-946.

## Examples

``` r
library(JWileymisc)
data(aces_daily, package = "JWileymisc")
m1 <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID),
  data = aces_daily)
m2 <- lme4::lmer(PosAff ~ STRESS + (1 | UserID),
  data = aces_daily)
m3 <- lme4::lmer(PosAff ~ STRESS + Female + (1 | UserID),
  data = aces_daily)

modelCompare(m1, m2)
#> $Comparison
#>         Model Estimator N_Obs     N_Groups        AIC        BIC         LL
#>        <char>    <char> <num>       <char>      <num>      <num>      <num>
#> 1:    Model 1      REML  6399 UserID (191) 13385.5617 13412.6173 -6688.7809
#> 2:    Model 2      REML  6399 UserID (191) 13184.7623 13225.3457 -6586.3812
#> 3: Difference               0               -200.7994  -187.2716   102.3997
#>     LLDF       Sigma  MarginalR2 ConditionalR2  MarginalF2 ConditionalF2
#>    <num>       <num>       <num>         <num>       <num>         <num>
#> 1:     4  0.65017805 0.142930809     0.6191536 0.166766944    1.62573070
#> 2:     6  0.63112686 0.148011774     0.6381639 0.173725140    1.76368250
#> 3:     2 -0.01905119 0.005080965     0.0190103 0.005963656    0.05253844
#>        Chi2            P
#>       <num>        <num>
#> 1:       NA           NA
#> 2:       NA           NA
#> 3: 204.7994 3.375832e-45
#> 
#> attr(,"class")
#> [1] "modelCompare.merMod" "modelCompare"       
modelCompare(m2, m3)
#> When using REML, the fixed effects structure must be identical,
#> but was different. Refitting with ML.
#> $Comparison
#>         Model Estimator N_Obs     N_Groups          AIC          BIC
#>        <char>    <char> <num>       <char>        <num>        <num>
#> 1:    Model 1        ML  6399 UserID (191) 13372.308546 13399.364134
#> 2:    Model 2        ML  6399 UserID (191) 13373.859321 13407.678806
#> 3: Difference               0                  1.550775     8.314672
#>               LL  LLDF         Sigma  MarginalR2 ConditionalR2  MarginalF2
#>            <num> <num>         <num>       <num>         <num>       <num>
#> 1: -6682.1542731     4  6.501264e-01 0.143317066  0.6182064065 0.167293009
#> 2: -6681.9296605     5  6.501263e-01 0.145209578  0.6185439889 0.169877404
#> 3:     0.2246126     1 -1.366298e-07 0.001892512  0.0003375823 0.002214008
#>    ConditionalF2      Chi2         P
#>            <num>     <num>     <num>
#> 1:  1.6192162915        NA        NA
#> 2:  1.6215342552        NA        NA
#> 3:  0.0008849837 0.4492251 0.5027032
#> 
#> attr(,"class")
#> [1] "modelCompare.merMod" "modelCompare"       

rm(m1, m2, m3)
```
