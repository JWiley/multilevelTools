# estimate detailed results per variable and effect sizes for both fixed and random effects from lmer models

This function extends the current
[`drop1`](https://rdrr.io/r/stats/add1.html) method for `merMod` class
objects from the lme4 package. Where the default method to be able to
drop both fixed and random effects at once.

## Usage

``` r
# S3 method for class 'merMod'
modelTest(object, method = c("Wald", "profile", "boot"), control, ...)
```

## Arguments

- object:

  A `link[lme4]{merMod-class}` object, the fitted result of
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

- method:

  A character vector indicating the types of confidence intervals to
  calculate. One of “Wald”, “profile”, or “boot”.

- control:

  A [`lmerControl`](https://rdrr.io/pkg/lme4/man/lmerControl.html)
  results used to control how models are estimated when updating.

- ...:

  Additional arguments passed to `confint`

## Details

At the moment, the function is aimed to
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) models and has very few
features for [`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) or
[`nlmer`](https://rdrr.io/pkg/lme4/man/nlmer.html) models. The primary
motivation was to provide a way to provide an overall test of whether a
variable “matters”. In multilevel data, a variable may be included in
both the fixed and random effects. To provide an overall test of whether
it matters requires jointly testing the fixed and random effects. This
also is needed to provide an overall effect size.

The function works by generating a formula with one specific variable or
“term” removed at all levels. A model is then fit on this reduced
formula and compared to the full model passed in. This is a complex
operation for mixed effects models for several reasons. Firstly, `R` has
no default mechanism for dropping terms from both the fixed and random
portions. Secondly, mixed effects models do not accomodate all types of
models. For example, if a model includes only a random slope with no
random intercept, if the random slope was dropped, there would be no
more random effects, and at that point,
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) or
[`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) will not run the
model. It is theoretically possible to instead fit the model using
[`lm`](https://rdrr.io/r/stats/lm.html) or
[`glm`](https://rdrr.io/r/stats/glm.html) but this becomes more complex
for certain model comparisons and calculations and is not currently
implemented. Marginal and conditional R2 values are calculated for each
term, and these are used also to calculate something akin to an
f-squared effect size.

This is a new function and it is important to carefully evaluate the
results and check that they are accurate and that they are sensible.
Check accuracy by viewing the model formulae for each reduced model and
checking that those are indeed accurate. In terms of checking whether a
result is sensible or not, there is a large literature on the difficulty
interpretting main effect tests in the presence of interactions. As it
is challenging to detect all interactions, especially ones that are made
outside of `R` formulae, all terms are tested. However, it likely does
not make sense to report results from dropping a main effect but keeping
the interaction term, so present and interpret these with caution.

## Examples

``` r
## these examples are slow to run
library(JWileymisc)
m1 <- lme4::lmer(extra ~ group + (1 | ID),
  data = sleep, REML = FALSE)
modelTest(m1)
#> $FixedEffects
#> Key: <Term>
#>           Term   Est   Type         LL       UL         Pval
#>         <char> <num> <char>      <num>    <num>        <num>
#> 1: (Intercept)  0.75     FE -0.3663711 1.866371 1.879250e-01
#> 2:      group2  1.58     FE  0.8567759 2.303224 1.853388e-05
#> 
#> $RandomEffects
#> Key: <Term>
#>                 Term       Est   Type    LL    UL
#>               <char>     <num> <char> <num> <num>
#> 1: sd_(Intercept)|ID 1.6010934     RE    NA    NA
#> 2:             sigma 0.8251061     RE    NA    NA
#> 
#> $EffectSizes
#>    Variable Estimator N_Obs N_Groups       AIC       BIC       LL  LLDF
#>      <char>    <char> <num>   <char>     <num>     <num>    <num> <num>
#> 1:    group               0          -8.414884 -7.419152 5.207442     1
#>         Sigma MarginalR2 ConditionalR2 MarginalF2 ConditionalF2     Chi2
#>         <num>      <num>         <num>      <num>         <num>    <num>
#> 1: -0.5637784  0.1613328     0.3226657  0.1923682      1.833431 10.41488
#>              P     FE     RE              Formula   Type
#>          <num> <lgcl> <lgcl>               <char> <fctr>
#> 1: 0.001250037   TRUE  FALSE extra ~ 1 + (1 | ID)  Fixed
#> 
#> $OverallModel
#> $Performance
#>     Model Estimator N_Obs N_Groups      AIC      BIC        LL  LLDF     Sigma
#>    <char>    <char> <num>   <char>    <num>    <num>     <num> <num>     <num>
#> 1: merMod        ML    20  ID (10) 78.50469 82.48762 -35.25235     4 0.8251061
#>    MarginalR2 ConditionalR2 MarginalF2 ConditionalF2
#>         <num>         <num>      <num>         <num>
#> 1:  0.1613328     0.8240099  0.1923682      4.682139
#> 
#> attr(,"class")
#> [1] "modelPerformance.merMod" "modelPerformance"       
#> 
#> attr(,"class")
#> [1] "modelTest.merMod" "modelTest"       

# \donttest{
data(aces_daily, package = "JWileymisc")

strictControl <- lme4::lmerControl(optCtrl = list(
   algorithm = "NLOPT_LN_NELDERMEAD",
   xtol_abs = 1e-10,
   ftol_abs = 1e-10))

m1 <- lme4::lmer(NegAff ~ STRESS + (1 + STRESS | UserID),
  data = aces_daily,
  control = strictControl)
modelTest(m1, method = "profile")
#> Computing profile confidence intervals ...
#> Parameters and CIs are based on REML, 
#> but modelTests requires ML not REML fit for comparisons, 
#> and these are used in effect sizes. Refitting.
#> $FixedEffects
#> Key: <Term>
#>           Term       Est   Type        LL        UL  Pval
#>         <char>     <num> <char>     <num>     <num> <num>
#> 1: (Intercept) 1.2357015     FE 1.2024389 1.2692167     0
#> 2:      STRESS 0.1183131     FE 0.1068207 0.1296297     0
#> 
#> $RandomEffects
#> Key: <Term>
#>                             Term        Est   Type         LL         UL
#>                           <char>      <num> <char>      <num>      <num>
#> 1: cor_STRESS.(Intercept)|UserID 0.47072409     RE 0.28594506 0.63920509
#> 2:         sd_(Intercept)|UserID 0.20585894     RE 0.17824126 0.23575940
#> 3:              sd_STRESS|UserID 0.06769381     RE 0.05902194 0.07720673
#> 4:                         sigma 0.41789606     RE 0.41053223 0.42549600
#> 
#> $EffectSizes
#>    Variable Estimator N_Obs N_Groups        AIC        BIC        LL  LLDF
#>      <char>    <char> <num>   <char>      <num>      <num>     <num> <num>
#> 1:   STRESS               0          -2701.2231 -2680.9361 1353.6115     3
#> 2:   STRESS               0           -538.1871  -524.6625  271.0936     2
#>          Sigma  MarginalR2 ConditionalR2  MarginalF2 ConditionalF2      Chi2
#>          <num>       <num>         <num>       <num>         <num>     <num>
#> 1: -0.10192041  0.22744094    0.11496664  0.29439942    0.25608128 2707.2231
#> 2: -0.02275066 -0.04844019    0.01044792 -0.06270096    0.02327211  542.1871
#>                P     FE     RE                            Formula
#>            <num> <lgcl> <lgcl>                             <char>
#> 1:  0.000000e+00   TRUE   TRUE          NegAff ~ 1 + (1 | UserID)
#> 2: 1.843162e-118  FALSE   TRUE NegAff ~ 1 + STRESS + (1 | UserID)
#>              Type
#>            <fctr>
#> 1: Fixed + Random
#> 2:         Random
#> 
#> $OverallModel
#> $Performance
#>     Model Estimator N_Obs     N_Groups      AIC      BIC        LL  LLDF
#>    <char>    <char> <num>       <char>    <num>    <num>     <num> <num>
#> 1: merMod        ML  6389 UserID (191) 7702.815 7743.389 -3845.407     6
#>        Sigma MarginalR2 ConditionalR2 MarginalF2 ConditionalF2
#>        <num>      <num>         <num>      <num>         <num>
#> 1: 0.4179028  0.2274409     0.5510541  0.2943994       1.22744
#> 
#> attr(,"class")
#> [1] "modelPerformance.merMod" "modelPerformance"       
#> 
#> attr(,"class")
#> [1] "modelTest.merMod" "modelTest"       

m2 <- lme4::lmer(NegAff ~ STRESS + I(STRESS^2) + (1 + STRESS | UserID),
  data = aces_daily, control = strictControl)

## might normally use more bootstraps but keeping low for faster run
modelTest(m2, method = "boot", nsim = 100)
#> Computing bootstrap confidence intervals ...
#> Error in eval(expr, p): object 'strictControl' not found
# }
```
