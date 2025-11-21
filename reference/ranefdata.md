# Create data and plots for [`brm`](https://paulbuerkner.com/brms/reference/brm.html) random effect models

Create data and plots for
[`brm`](https://paulbuerkner.com/brms/reference/brm.html) random effect
models

## Usage

``` r
ranefdata(object, usevars, newdata, idvar, CI = 0.95, robust = FALSE)
```

## Arguments

- object:

  a
  [`brmsfit-class`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object

- usevars:

  a character vector of random effects to plot

- newdata:

  a `data.table` object with the data used to generate the random
  effects, this is used as an anchor for the random intercepts so they
  have a meaningful 0 point

- idvar:

  a character string specifying the grouping variable name for the
  random effects

- CI:

  a numeric value between 0 and 1 specifying the interval to use.
  Defaults to 0.95.

- robust:

  a logical value indicating whether to use robust estimates or not.
  Defaults to FALSE. Passed on to
  [`posterior_summary`](https://paulbuerkner.com/brms/reference/posterior_summary.html)
  and
  [`.summary.ID`](https://joshuawiley.com/multilevelTools/reference/dot-summary.ID.md).

## Value

a list with the following components: \* `plot`: a list of ggplot
objects \* `plotdat`: a list of data.table objects with the data used to
generate the plots \* `relong`: a data.table object with the random
effects in long format \* `yhat`: a list of data.table objects with the
expected values for the random effects \* `usevars`: a character vector
of the random effects to plot \* `idvar`: a character string specifying
the grouping variable name for the random effects

## Examples

``` r
# \donttest{
if (FALSE) {
library(data.table)
library(brms)
library(ggpubr)

current.seed <- .Random.seed
set.seed(12345)
nGroups <- 100
nObs <- 20
theta.location <- matrix(rnorm(nGroups * 2), nrow = nGroups, ncol = 2)
theta.location[, 1] <- theta.location[, 1] - mean(theta.location[, 1])
theta.location[, 2] <- theta.location[, 2] - mean(theta.location[, 2])
theta.location[, 1] <- theta.location[, 1] / sd(theta.location[, 1])
theta.location[, 2] <- theta.location[, 2] / sd(theta.location[, 2])
theta.location <- theta.location %*% chol(matrix(c(1.5, -.25, -.25, .5^2), 2))
theta.location[, 1] <- theta.location[, 1] - 2.5
theta.location[, 2] <- theta.location[, 2] + 1

dmixed <- data.table(
  x = rep(rep(0:1, each = nObs / 2), times = nGroups))
  dmixed[, ID := rep(seq_len(nGroups), each = nObs)]

  for (i in seq_len(nGroups)) {
    dmixed[ID == i, y := rnorm(
      n = nObs,
      mean = theta.location[i, 1] + theta.location[i, 2] * x,
      sd = exp(1 + theta.location[i, 1] + theta.location[i, 2] * x))
      ]
  }

## note this model takes several minutes, even on a high performance machine
ls.me <- brm(bf(
  y ~ 1 + x + (1 + x | p | ID),
  sigma ~ 1 + x + (1 + x | p | ID)),
  family = "gaussian",
  data = dmixed, seed = 1234,
  silent = 2, refresh = 0, iter = 2000, warmup = 1000, thin = 1,
  chains = 4L, cores = 4L)

out <- ranefdata(
  ls.me,
  newdata = data.table(ID = dmixed$ID[1], x = 0), 
  usevars = c("Intercept", "x", "sigma_Intercept", "sigma_x"),
  idvar = "ID")

do.call(ggarrange, c(out$replots, ncol=2,nrow=2))
do.call(ggarrange, c(out$scatterplots, ncol=2,nrow=3))

## set seed back to what it was
set.seed(current.seed)

## cleanup
rm(current.seed, nGroups, nObs, theta.location, dmixed, ls.me, out)
}
# }
```
