# Format results from a linear mixed model

Format results from a linear mixed model

## Usage

``` r
# S3 method for class 'modelTest.merMod'
APAStyler(
  object,
  format = list(FixedEffects = c("%s%s [%s, %s]"), RandomEffects = c("%s",
    "%s [%s, %s]"), EffectSizes = c("%s/%s, %s")),
  digits = 2,
  pcontrol = list(digits = 3, stars = TRUE, includeP = FALSE, includeSign = FALSE,
    dropLeadingZero = TRUE),
  ...
)
```

## Arguments

- object:

  A list of one (or more) models estimated from lmer

- format:

  A list giving the formatting style to be used for the fixed effects,
  random effects, and effect sizes. For the random effects, must be two
  options, one for when the random effects do not have confidence
  intervals and one when the random effects do have confidence
  intervals.

- digits:

  A numeric value indicating the number of digits to print. This is
  still in early implementation stages and currently does not change all
  parts of the output (which default to 2 decimals per APA style).

- pcontrol:

  A list controlling how p values are formatted.

- ...:

  Additional arguments passed to
  [`confint`](https://rdrr.io/r/stats/confint.html). Notably `nsim` and
  `boot.type` if the bootstrap method is used.

## Value

a data table of character data

## Examples

``` r
library(JWileymisc)
data(sleepstudy, package = "lme4")

m1 <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
  data = sleepstudy)
m2 <- lme4::lmer(Reaction ~ Days + I(Days^2) + (1 + Days | Subject),
  data = sleepstudy)

testm1 <- modelTest(m1)
#> Parameters and CIs are based on REML, 
#> but modelTests requires ML not REML fit for comparisons, 
#> and these are used in effect sizes. Refitting.
testm2 <- modelTest(m2)
#> Parameters and CIs are based on REML, 
#> but modelTests requires ML not REML fit for comparisons, 
#> and these are used in effect sizes. Refitting.

APAStyler(testm1)
#>                             Term                        Est           Type
#>                           <char>                     <char>         <char>
#>  1:                  (Intercept) 251.41*** [238.03, 264.78]  Fixed Effects
#>  2:                         Days  10.47*** [  7.44,  13.50]  Fixed Effects
#>  3: cor_Days.(Intercept)|Subject                       0.07 Random Effects
#>  4:       sd_(Intercept)|Subject                      24.74 Random Effects
#>  5:              sd_Days|Subject                       5.92 Random Effects
#>  6:                        sigma                      25.59 Random Effects
#>  7:                     Model DF                          6  Overall Model
#>  8:                   N (Groups)               Subject (18)  Overall Model
#>  9:             N (Observations)                        180  Overall Model
#> 10:                       logLik                    -875.97  Overall Model
#> 11:                          AIC                    1763.94  Overall Model
#> 12:                          BIC                    1783.10  Overall Model
#> 13:                  Marginal R2                       0.29  Overall Model
#> 14:                  Marginal F2                       0.40  Overall Model
#> 15:               Conditional R2                       0.79  Overall Model
#> 16:               Conditional F2                       3.82  Overall Model
#> 17:        Days (Fixed + Random)        0.40/1.99, p < .001   Effect Sizes
#> 18:                Days (Random)        0.00/0.46, p < .001   Effect Sizes
APAStyler(list(Linear = testm1, Quadratic = testm2))
#>                             Term                     Linear
#>                           <char>                     <char>
#>  1:                  (Intercept) 251.41*** [238.03, 264.78]
#>  2:                         Days  10.47*** [  7.44,  13.50]
#>  3:                    I(Days^2)                           
#>  4: cor_Days.(Intercept)|Subject                       0.07
#>  5:       sd_(Intercept)|Subject                      24.74
#>  6:              sd_Days|Subject                       5.92
#>  7:                        sigma                      25.59
#>  8:                     Model DF                          6
#>  9:                   N (Groups)               Subject (18)
#> 10:             N (Observations)                        180
#> 11:                       logLik                    -875.97
#> 12:                          AIC                    1763.94
#> 13:                          BIC                    1783.10
#> 14:                  Marginal R2                       0.29
#> 15:                  Marginal F2                       0.40
#> 16:               Conditional R2                       0.79
#> 17:               Conditional F2                       3.82
#> 18:        Days (Fixed + Random)        0.40/1.99, p < .001
#> 19:                Days (Random)        0.00/0.46, p < .001
#> 20:            I(Days^2) (Fixed)                           
#>                             Term                     Linear
#>                      Quadratic           Type
#>                         <char>         <char>
#>  1: 255.45*** [240.72, 270.18]  Fixed Effects
#>  2:    7.43** [  1.91,  12.96]  Fixed Effects
#>  3:      0.34 [ -0.18,   0.85]  Fixed Effects
#>  4:                       0.06 Random Effects
#>  5:                      24.76 Random Effects
#>  6:                       5.93 Random Effects
#>  7:                      25.53 Random Effects
#>  8:                          7  Overall Model
#>  9:               Subject (18)  Overall Model
#> 10:                        180  Overall Model
#> 11:                    -875.14  Overall Model
#> 12:                    1764.28  Overall Model
#> 13:                    1786.63  Overall Model
#> 14:                       0.29  Overall Model
#> 15:                       0.41  Overall Model
#> 16:                       0.79  Overall Model
#> 17:                       3.87  Overall Model
#> 18:        0.01/0.52, p < .001   Effect Sizes
#> 19:        0.00/0.46, p < .001   Effect Sizes
#> 20:        0.00/0.01, p = .198   Effect Sizes
#>                      Quadratic           Type
APAStyler(testm1,
  format = list(
    FixedEffects = "%s, %s (%s, %s)",
    RandomEffects = c("%s", "%s (%s, %s)"),
    EffectSizes = "%s, %s; %s"),
  pcontrol = list(digits = 3, stars = FALSE,
                  includeP = TRUE, includeSign = TRUE,
                  dropLeadingZero = TRUE))
#>                             Term                               Est
#>                           <char>                            <char>
#>  1:                  (Intercept) 251.41, p < .001 (238.03, 264.78)
#>  2:                         Days  10.47, p < .001 (  7.44,  13.50)
#>  3: cor_Days.(Intercept)|Subject                              0.07
#>  4:       sd_(Intercept)|Subject                             24.74
#>  5:              sd_Days|Subject                              5.92
#>  6:                        sigma                             25.59
#>  7:                     Model DF                                 6
#>  8:                   N (Groups)                      Subject (18)
#>  9:             N (Observations)                               180
#> 10:                       logLik                           -875.97
#> 11:                          AIC                           1763.94
#> 12:                          BIC                           1783.10
#> 13:                  Marginal R2                              0.29
#> 14:                  Marginal F2                              0.40
#> 15:               Conditional R2                              0.79
#> 16:               Conditional F2                              3.82
#> 17:        Days (Fixed + Random)              0.40, 1.99; p < .001
#> 18:                Days (Random)              0.00, 0.46; p < .001
#>               Type
#>             <char>
#>  1:  Fixed Effects
#>  2:  Fixed Effects
#>  3: Random Effects
#>  4: Random Effects
#>  5: Random Effects
#>  6: Random Effects
#>  7:  Overall Model
#>  8:  Overall Model
#>  9:  Overall Model
#> 10:  Overall Model
#> 11:  Overall Model
#> 12:  Overall Model
#> 13:  Overall Model
#> 14:  Overall Model
#> 15:  Overall Model
#> 16:  Overall Model
#> 17:   Effect Sizes
#> 18:   Effect Sizes

# \donttest{

testm1 <- modelTest(m1, method = "profile")
#> Computing profile confidence intervals ...
#> Parameters and CIs are based on REML, 
#> but modelTests requires ML not REML fit for comparisons, 
#> and these are used in effect sizes. Refitting.
testm2 <- modelTest(m2, method = "profile")
#> Computing profile confidence intervals ...
#> Parameters and CIs are based on REML, 
#> but modelTests requires ML not REML fit for comparisons, 
#> and these are used in effect sizes. Refitting.

APAStyler(testm1)
#>                             Term                        Est           Type
#>                           <char>                     <char>         <char>
#>  1:                  (Intercept) 251.41*** [237.68, 265.13]  Fixed Effects
#>  2:                         Days  10.47*** [  7.36,  13.58]  Fixed Effects
#>  3: cor_Days.(Intercept)|Subject        0.07 [-0.48,  0.68] Random Effects
#>  4:       sd_(Intercept)|Subject       24.74 [14.38, 37.72] Random Effects
#>  5:              sd_Days|Subject        5.92 [ 3.80,  8.75] Random Effects
#>  6:                        sigma       25.59 [22.90, 28.86] Random Effects
#>  7:                     Model DF                          6  Overall Model
#>  8:                   N (Groups)               Subject (18)  Overall Model
#>  9:             N (Observations)                        180  Overall Model
#> 10:                       logLik                    -875.97  Overall Model
#> 11:                          AIC                    1763.94  Overall Model
#> 12:                          BIC                    1783.10  Overall Model
#> 13:                  Marginal R2                       0.29  Overall Model
#> 14:                  Marginal F2                       0.40  Overall Model
#> 15:               Conditional R2                       0.79  Overall Model
#> 16:               Conditional F2                       3.82  Overall Model
#> 17:        Days (Fixed + Random)        0.40/1.99, p < .001   Effect Sizes
#> 18:                Days (Random)        0.00/0.46, p < .001   Effect Sizes
APAStyler(list(Linear = testm1, Quadratic = testm2))
#>                             Term                     Linear
#>                           <char>                     <char>
#>  1:                  (Intercept) 251.41*** [237.68, 265.13]
#>  2:                         Days  10.47*** [  7.36,  13.58]
#>  3:                    I(Days^2)                           
#>  4: cor_Days.(Intercept)|Subject        0.07 [-0.48,  0.68]
#>  5:       sd_(Intercept)|Subject       24.74 [14.38, 37.72]
#>  6:              sd_Days|Subject        5.92 [ 3.80,  8.75]
#>  7:                        sigma       25.59 [22.90, 28.86]
#>  8:                     Model DF                          6
#>  9:                   N (Groups)               Subject (18)
#> 10:             N (Observations)                        180
#> 11:                       logLik                    -875.97
#> 12:                          AIC                    1763.94
#> 13:                          BIC                    1783.10
#> 14:                  Marginal R2                       0.29
#> 15:                  Marginal F2                       0.40
#> 16:               Conditional R2                       0.79
#> 17:               Conditional F2                       3.82
#> 18:        Days (Fixed + Random)        0.40/1.99, p < .001
#> 19:                Days (Random)        0.00/0.46, p < .001
#> 20:            I(Days^2) (Fixed)                           
#>                             Term                     Linear
#>                      Quadratic           Type
#>                         <char>         <char>
#>  1: 255.45*** [240.52, 270.38]  Fixed Effects
#>  2:    7.43** [  1.93,  12.94]  Fixed Effects
#>  3:      0.34 [ -0.18,   0.85]  Fixed Effects
#>  4:        0.06 [-0.48,  0.68] Random Effects
#>  5:       24.76 [14.48, 37.75] Random Effects
#>  6:        5.93 [ 3.81,  8.76] Random Effects
#>  7:       25.53 [22.77, 28.69] Random Effects
#>  8:                          7  Overall Model
#>  9:               Subject (18)  Overall Model
#> 10:                        180  Overall Model
#> 11:                    -875.14  Overall Model
#> 12:                    1764.28  Overall Model
#> 13:                    1786.63  Overall Model
#> 14:                       0.29  Overall Model
#> 15:                       0.41  Overall Model
#> 16:                       0.79  Overall Model
#> 17:                       3.87  Overall Model
#> 18:        0.01/0.52, p < .001   Effect Sizes
#> 19:        0.00/0.46, p < .001   Effect Sizes
#> 20:        0.00/0.01, p = .198   Effect Sizes
#>                      Quadratic           Type
APAStyler(testm1,
  format = list(
    FixedEffects = "%s, %s (%s, %s)",
    RandomEffects = c("%s", "%s (%s, %s)"),
    EffectSizes = "%s, %s; %s"),
  pcontrol = list(digits = 3, stars = FALSE,
                  includeP = TRUE, includeSign = TRUE,
                  dropLeadingZero = TRUE))
#>                             Term                               Est
#>                           <char>                            <char>
#>  1:                  (Intercept) 251.41, p < .001 (237.68, 265.13)
#>  2:                         Days  10.47, p < .001 (  7.36,  13.58)
#>  3: cor_Days.(Intercept)|Subject               0.07 (-0.48,  0.68)
#>  4:       sd_(Intercept)|Subject              24.74 (14.38, 37.72)
#>  5:              sd_Days|Subject               5.92 ( 3.80,  8.75)
#>  6:                        sigma              25.59 (22.90, 28.86)
#>  7:                     Model DF                                 6
#>  8:                   N (Groups)                      Subject (18)
#>  9:             N (Observations)                               180
#> 10:                       logLik                           -875.97
#> 11:                          AIC                           1763.94
#> 12:                          BIC                           1783.10
#> 13:                  Marginal R2                              0.29
#> 14:                  Marginal F2                              0.40
#> 15:               Conditional R2                              0.79
#> 16:               Conditional F2                              3.82
#> 17:        Days (Fixed + Random)              0.40, 1.99; p < .001
#> 18:                Days (Random)              0.00, 0.46; p < .001
#>               Type
#>             <char>
#>  1:  Fixed Effects
#>  2:  Fixed Effects
#>  3: Random Effects
#>  4: Random Effects
#>  5: Random Effects
#>  6: Random Effects
#>  7:  Overall Model
#>  8:  Overall Model
#>  9:  Overall Model
#> 10:  Overall Model
#> 11:  Overall Model
#> 12:  Overall Model
#> 13:  Overall Model
#> 14:  Overall Model
#> 15:  Overall Model
#> 16:  Overall Model
#> 17:   Effect Sizes
#> 18:   Effect Sizes

# }

rm(m1, m2, testm1, testm2)
```
