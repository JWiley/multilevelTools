# modelPerformance method for merMod objects

For pseudo R2 by cluster, the squared correlation between observed and
predicted values for each cluster unit is returned. For the overall
model, the marginal and conditional R2 are calculated as described in
the references.

## Usage

``` r
# S3 method for class 'merMod'
modelPerformance(object, ...)
```

## Arguments

- object:

  A model from `lmer`.

- ...:

  Additional arguments, not currently used.

## Value

a named vector with the marginal and conditional R2 values, if
`CLUSTER = FALSE`, otherwise, a data table with the pseudo R2 for each
cluster unit. A list with a `data.table` with the following elements: -
Model: A character string indicating the model type, here merMod -
Estimator: A character string indicating whether the model was estimated
with REML or ML - N_Obs: The number of observations - N_Groups: A
character string indicating the number of unique units in each
grouping/clustering variable. - AIC: Akaike Information Criterion - BIC:
Bayesian Information Criterion - LL: log likelihood - LLDF: log
likelihood degrees of freedom - Sigma: Residual standard deviation -
MarginalR2: in sample variance explained by the fixed effects -
ConditionalR2: in sample variance explained by the fixed and random
effects - MarginalF2: Cohen's F2 effect size R2 / (1 - R2) based off the
Marginal R2 - ConditionalF2: Cohen's F2 effect size R2 / (1 - R2) based
off the Conditional R2

## References

For estimating the marginal and conditional R-squared values, see:
Nakagawa, S. and Schielzeth, H. (2013)
\<doi:10.1111/j.2041-210x.2012.00261.x\> "A general and simple method
for obtaining R2 from generalized linear mixed-effects models" and also:
Johnson, P. C. (2014) \<doi:10.1111/2041-210X.12225\> "Extension of
Nakagawa & Schielzeth's R2GLMM to random slopes models"

## Examples

``` r
library(JWileymisc)
data(aces_daily, package = "JWileymisc")
m1 <- lme4::lmer(PosAff ~ 1 + (1 | UserID),
  data = aces_daily)
modelPerformance(m1)
#> $Performance
#>     Model Estimator N_Obs     N_Groups      AIC      BIC        LL  LLDF
#>    <char>    <char> <num>       <char>    <num>    <num>     <num> <num>
#> 1: merMod      REML  6399 UserID (191) 14800.56 14820.85 -7397.281     3
#>        Sigma MarginalR2 ConditionalR2 MarginalF2 ConditionalF2
#>        <num>      <num>         <num>      <num>         <num>
#> 1: 0.7273136          0     0.5431473          0       1.18889
#> 
#> attr(,"class")
#> [1] "modelPerformance.merMod" "modelPerformance"       

m1 <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID),
  data = aces_daily)
modelPerformance(m1)
#> $Performance
#>     Model Estimator N_Obs     N_Groups      AIC      BIC        LL  LLDF
#>    <char>    <char> <num>       <char>    <num>    <num>     <num> <num>
#> 1: merMod      REML  6399 UserID (191) 13184.76 13225.35 -6586.381     6
#>        Sigma MarginalR2 ConditionalR2 MarginalF2 ConditionalF2
#>        <num>      <num>         <num>      <num>         <num>
#> 1: 0.6311269  0.1480118     0.6381639  0.1737251      1.763683
#> 
#> attr(,"class")
#> [1] "modelPerformance.merMod" "modelPerformance"       

rm(m1)
```
