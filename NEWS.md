# multilevelTools 0.1.3

## Bug Fixes

* `modelTest()` no longer fails for models with a continuous x categorical interaction.
  Estimates for dropping the "simple" effect of the continuous variable are still 
  not calculable, but the rest of the calculations are still performed and that line
  is simply set to NA.
  
## Changes

* moved to testthat 3rd edition
* moved CI to GitHub actions
* use preferably for package website

# multilevelTools 0.1.2

## New Features

* New `weighted.sma` function to calculate weighted simple 
  moving averages.

# multilevelTools 0.1.1

## New Features

* Beta methods to support lme models, class `lme` for 
  `residualDiagnostics()` and `modelDiagnostics()`
  with more planned in future updates.

# multilevelTools 0.1.0

## New Features

* Methods to support lme4 models, class `merMod` for 
  `modelTest()`, `modelDiagnostics()`, and `APAStyler()`.

* New vignette added showing sample use case of the package.

## Ported Features

* `omegaSEM()` Function that calculates coefficient omega for 
  measuring internal consistency reliability. Works for two 
  level models and returns within and between level omega 
  values.
  
* `R2.merMod()` A method to calculate the marginal and
  conditional variance accounted for by a model estimated by 
  `lmer()`.

* `modelCompare.merMod()` A method to compare two models estimated by 
   `lmer()` include significance tests and effect sizes 
   for estimates of the variance explained.
          
* `iccMixed()` A function to calculate the intraclass correlation 
	coefficient using mixed effects models.  Works with either 
	normally distributed outcomes or binary outcomes, in which case 
	the latent variable estimate of the ICC is computed.
	
* `nEffective()` Calculates the effective sample size based on 
    the number of independent units, number of observations per 
	unit, and the intraclass correlation coefficient.
	
* `acfByID()` Calculates the lagged autocorrelation of a variable 
    by an ID variable and returns a data.table for further use,
	such as examination, summary, or plotting

* `meanDecompose()` function added to decompose multilevel or 
    repeated measures data into means and residuals.

* `meanDeviations()` A simple function to calculate means and mean 
	deviations, useful for creating between and within versions of 
	a variable in a data.table
		
