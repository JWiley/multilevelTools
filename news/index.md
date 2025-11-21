# Changelog

## multilevelTools 0.2.1

CRAN release: 2025-09-07

### Changes

- Improved documentation links to functions in other packages.

## multilevelTools 0.2.0

CRAN release: 2025-04-13

### New Features

- [`evaluateLags()`](https://joshuawiley.com/multilevelTools/reference/evaluateLags.md)
  function to evaluate how far back lags should go to improve model fit.
- [`ranefdata()`](https://joshuawiley.com/multilevelTools/reference/ranefdata.md)
  function to extract random effects from a `brms` model and return and
  plot them. Designed to make caterpilar plots with posterior summaries.

### Changes

- Moved testthat to 3rd edition.
- Updated GitHub actions

## multilevelTools 0.1.3

### Bug Fixes

- [`modelTest()`](https://joshuawiley.com/JWileymisc/reference/modelTest.html)
  no longer fails for models with a continuous x categorical
  interaction. Estimates for dropping the “simple” effect of the
  continuous variable are still not calculable, but the rest of the
  calculations are still performed and that line is simply set to NA.

### Changes

- moved to testthat 3rd edition
- moved CI to GitHub actions
- use preferably for package website

## multilevelTools 0.1.2

### New Features

- New `weighted.sma` function to calculate weighted simple moving
  averages.

## multilevelTools 0.1.1

CRAN release: 2020-03-04

### New Features

- Beta methods to support lme models, class `lme` for
  [`residualDiagnostics()`](https://joshuawiley.com/JWileymisc/reference/residualDiagnostics.html)
  and
  [`modelDiagnostics()`](https://joshuawiley.com/JWileymisc/reference/modelDiagnostics.html)
  with more planned in future updates.

## multilevelTools 0.1.0

### New Features

- Methods to support lme4 models, class `merMod` for
  [`modelTest()`](https://joshuawiley.com/JWileymisc/reference/modelTest.html),
  [`modelDiagnostics()`](https://joshuawiley.com/JWileymisc/reference/modelDiagnostics.html),
  and
  [`APAStyler()`](https://joshuawiley.com/JWileymisc/reference/APAStyler.html).

- New vignette added showing sample use case of the package.

### Ported Features

- [`omegaSEM()`](https://joshuawiley.com/multilevelTools/reference/omegaSEM.md)
  Function that calculates coefficient omega for measuring internal
  consistency reliability. Works for two level models and returns within
  and between level omega values.

- [`R2.merMod()`](https://joshuawiley.com/multilevelTools/reference/R2.merMod.md)
  A method to calculate the marginal and conditional variance accounted
  for by a model estimated by
  [`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).

- [`modelCompare.merMod()`](https://joshuawiley.com/multilevelTools/reference/modelCompare.merMod.md)
  A method to compare two models estimated by
  [`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) include
  significance tests and effect sizes for estimates of the variance
  explained.

- [`iccMixed()`](https://joshuawiley.com/multilevelTools/reference/iccMixed.md)
  A function to calculate the intraclass correlation coefficient using
  mixed effects models. Works with either normally distributed outcomes
  or binary outcomes, in which case the latent variable estimate of the
  ICC is computed.

- [`nEffective()`](https://joshuawiley.com/multilevelTools/reference/nEffective.md)
  Calculates the effective sample size based on the number of
  independent units, number of observations per unit, and the intraclass
  correlation coefficient.

- [`acfByID()`](https://joshuawiley.com/multilevelTools/reference/acfByID.md)
  Calculates the lagged autocorrelation of a variable by an ID variable
  and returns a data.table for further use, such as examination,
  summary, or plotting

- [`meanDecompose()`](https://joshuawiley.com/multilevelTools/reference/meanDecompose.md)
  function added to decompose multilevel or repeated measures data into
  means and residuals.

- [`meanDeviations()`](https://joshuawiley.com/multilevelTools/reference/meanDeviations.md)
  A simple function to calculate means and mean deviations, useful for
  creating between and within versions of a variable in a data.table
