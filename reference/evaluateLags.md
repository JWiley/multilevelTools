# Create lag variables and evaluate models with different number of lags

This function creates the desired number of lags and tests consecutive
models from a model with no lags (lag 0), lag 0 + lag1, etc. and reports
model performance. This helps evaluate how many lags are needed.

## Usage

``` r
evaluateLags(formula, lagvar, nlags = 0L, idvar, data, ...)
```

## Arguments

- formula:

  A `character` string giving the
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) formula to use as a
  base. The variable to be tested with lags gets added as fixed effects
  only to this, currently.

- lagvar:

  A `character` string giving the name of the variable to test lags for.

- nlags:

  An `integer` (e.g., 0L, 3L) giving the number of lags to test.
  Defaults to 0L but really should be more. Must be a positive integer.

- idvar:

  A `character` string giving the name o the ID variable.

- data:

  A `data.table` dataset ideally or at least a `data.frame`.

- ...:

  Additional arguments passed to `lmer`, used to control model fitting.

## Details

Currently only linear mixed effects models are allowed.

## Examples

``` r
## these examples are slow to run
data(aces_daily, package = "JWileymisc")

evaluateLags(
 "NegAff ~ Female + Age + BornAUS + (1 | UserID)",
 "STRESS",
 4L,
 "UserID",
 aces_daily)
#>      Lag  Model Estimator N_Obs     N_Groups      AIC      BIC        LL  LLDF
#>    <int> <char>    <char> <num>       <char>    <num>    <num>     <num> <num>
#> 1:     0 merMod        ML  5372 UserID (190) 6751.190 6797.313 -3368.595     7
#> 2:     1 merMod        ML  5372 UserID (190) 6681.458 6734.169 -3332.729     8
#> 3:     2 merMod        ML  5372 UserID (190) 6648.660 6707.960 -3315.330     9
#> 4:     3 merMod        ML  5372 UserID (190) 6638.343 6704.233 -3309.172    10
#> 5:     4 merMod        ML  5372 UserID (190) 6640.292 6712.771 -3309.146    11
#>        Sigma MarginalR2 ConditionalR2 MarginalF2 ConditionalF2
#>        <num>      <num>         <num>      <num>         <num>
#> 1: 0.4301553  0.2941101     0.5707243  0.4166516      1.329505
#> 2: 0.4276613  0.3233266     0.5801990  0.4778178      1.382081
#> 3: 0.4264497  0.3399432     0.5867621  0.5150211      1.419913
#> 4: 0.4260388  0.3483062     0.5903241  0.5344629      1.440954
#> 5: 0.4260312  0.3478455     0.5901544  0.5333790      1.439943
#>                                                                                                        Formula
#>                                                                                                         <char>
#> 1:                                                     NegAff ~ Female + Age + BornAUS + (1 | UserID) + STRESS
#> 2:                                        NegAff ~ Female + Age + BornAUS + (1 | UserID) + STRESS + STRESSLag1
#> 3:                           NegAff ~ Female + Age + BornAUS + (1 | UserID) + STRESS + STRESSLag1 + STRESSLag2
#> 4:              NegAff ~ Female + Age + BornAUS + (1 | UserID) + STRESS + STRESSLag1 + STRESSLag2 + STRESSLag3
#> 5: NegAff ~ Female + Age + BornAUS + (1 | UserID) + STRESS + STRESSLag1 + STRESSLag2 + STRESSLag3 + STRESSLag4

# \donttest{
 ## not run, more complex example with random slope, fails to converge
evaluateLags(
 "NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID)",
 "STRESS",
 5L,
 "UserID",
 aces_daily)
#> Warning: Model failed to converge with max|grad| = 0.00302696 (tol = 0.002, component 1)
#>      Lag  Model Estimator N_Obs     N_Groups      AIC      BIC        LL  LLDF
#>    <int> <char>    <char> <num>       <char>    <num>    <num>     <num> <num>
#> 1:     0 merMod        ML  5144 UserID (190) 5903.510 5962.420 -2942.755     9
#> 2:     1 merMod        ML  5144 UserID (190) 5837.591 5903.046 -2908.795    10
#> 3:     2 merMod        ML  5144 UserID (190) 5816.798 5888.799 -2897.399    11
#> 4:     3 merMod        ML  5144 UserID (190) 5806.522 5885.069 -2891.261    12
#> 5:     4 merMod        ML  5144 UserID (190) 5808.423 5893.515 -2891.211    13
#> 6:     5 merMod        ML  5144 UserID (190) 5810.158 5901.796 -2891.079    14
#>        Sigma MarginalR2 ConditionalR2 MarginalF2 ConditionalF2
#>        <num>      <num>         <num>      <num>         <num>
#> 1: 0.4021802  0.2441862     0.5872656  0.3230772      1.422865
#> 2: 0.4001612  0.2785727     0.5956186  0.3861410      1.472913
#> 3: 0.3994679  0.2943976     0.6000910  0.4172287      1.500569
#> 4: 0.3990020  0.3040614     0.6039227  0.4369085      1.524760
#> 5: 0.3989923  0.3033384     0.6036954  0.4354171      1.523311
#> 6: 0.3989851  0.3044128     0.6041170  0.4376342      1.525999
#>                                                                                                                              Formula
#>                                                                                                                               <char>
#> 1:                                                                  NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS
#> 2:                                                     NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1
#> 3:                                        NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1 + STRESSLag2
#> 4:                           NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1 + STRESSLag2 + STRESSLag3
#> 5:              NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1 + STRESSLag2 + STRESSLag3 + STRESSLag4
#> 6: NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1 + STRESSLag2 + STRESSLag3 + STRESSLag4 + STRESSLag5

## use different control to fit model and now converges
strictControl <- lme4::lmerControl(optCtrl = list(
   algorithm = "NLOPT_LN_NELDERMEAD",
   xtol_abs = 1e-10,
   ftol_abs = 1e-10))
evaluateLags(
 "NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID)",
 "STRESS",
 5L,
 "UserID",
 aces_daily,
control = strictControl)
#>      Lag  Model Estimator N_Obs     N_Groups      AIC      BIC        LL  LLDF
#>    <int> <char>    <char> <num>       <char>    <num>    <num>     <num> <num>
#> 1:     0 merMod        ML  5144 UserID (190) 5903.510 5962.420 -2942.755     9
#> 2:     1 merMod        ML  5144 UserID (190) 5837.591 5903.046 -2908.795    10
#> 3:     2 merMod        ML  5144 UserID (190) 5816.798 5888.799 -2897.399    11
#> 4:     3 merMod        ML  5144 UserID (190) 5806.522 5885.069 -2891.261    12
#> 5:     4 merMod        ML  5144 UserID (190) 5808.423 5893.515 -2891.211    13
#> 6:     5 merMod        ML  5144 UserID (190) 5810.158 5901.796 -2891.079    14
#>        Sigma MarginalR2 ConditionalR2 MarginalF2 ConditionalF2
#>        <num>      <num>         <num>      <num>         <num>
#> 1: 0.4021799  0.2441830     0.5872689  0.3230717      1.422885
#> 2: 0.4001612  0.2785749     0.5956152  0.3861453      1.472892
#> 3: 0.3994694  0.2944043     0.6000779  0.4172423      1.500487
#> 4: 0.3990019  0.3040583     0.6039258  0.4369019      1.524779
#> 5: 0.3989926  0.3033362     0.6036964  0.4354127      1.523318
#> 6: 0.3989852  0.3044130     0.6041160  0.4376347      1.525992
#>                                                                                                                              Formula
#>                                                                                                                               <char>
#> 1:                                                                  NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS
#> 2:                                                     NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1
#> 3:                                        NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1 + STRESSLag2
#> 4:                           NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1 + STRESSLag2 + STRESSLag3
#> 5:              NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1 + STRESSLag2 + STRESSLag3 + STRESSLag4
#> 6: NegAff ~ Female + Age + BornAUS + (1 + STRESS | UserID) + STRESS + STRESSLag1 + STRESSLag2 + STRESSLag3 + STRESSLag4 + STRESSLag5
# }
```
