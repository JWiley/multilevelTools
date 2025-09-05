#' Mean decomposition of a variable by group(s)
#'
#' This function decomposes a variable in a long data set by grouping
#' factors, such as by ID. If multiple grouping factors are listed,
#' the decomposition is in order from left to right.
#' Residuals from the lowest level are returned.
#'
#' @param formula A formula of the variables to be used in the analysis.
#'   Should have the form: variable ~ groupingfactors.
#' @param data A data table or data frame containing the variables
#'   used in the formula.  This is a required argument.
#' @return A list of data tables with the means or residuals
#' @keywords multivariate
#' @importFrom stats terms
#' @importFrom data.table is.data.table as.data.table :=
#' @export
#' @examples
#' meanDecompose(mpg ~ vs, data = mtcars)
#' meanDecompose(mpg ~ vs + cyl, data = mtcars)
#'
#' ## Example plotting the results
#' tmp <- meanDecompose(Sepal.Length ~ Species, data = iris)
#' do.call(ggpubr::ggarrange, c(lapply(names(tmp), function(x) {
#'   plot(JWileymisc::testDistribution(tmp[[x]]$X), plot = FALSE, varlab = x)$Density
#' }), ncol = 1))
#'
#' rm(tmp)
meanDecompose <- function(formula, data) {
  v <- all.vars(formula)

  if (!is.data.table(data)) {
    data <- as.data.table(data)[, v, with = FALSE]
  } else {
    data <- data[, v, with = FALSE]
  }

  out <- vector("list", length = length(v))

  vres <- paste0(v[1], "_residual")
  stopifnot(!any(vres %in% v))

  data[, (vres) := get(v[1])]

  vfinal <- vector("character", length = length(v))

  for (i in 2:length(v)) {
    vname <- paste0(v[1], "_", v[i])
    data[, (vname) := mean(get(vres), na.rm = TRUE), by = c(v[2:i])]
    data[, (vres) := get(vres) - get(vname)]
    out[[i - 1]] <- data[, .(X = get(vname)[1]), by = c(v[2:i])]
    vfinal[i - 1] <- paste0(v[1], " by ", paste(v[2:i], collapse = " & "))
  }
  out[[length(v)]] <- data[, .(X = get(vres))]
  vfinal[length(v)] <- paste0(v[1], " by ", "residual")

  names(out) <- vfinal

  return(out)
}

# clear R CMD CHECK notes
if(getRversion() >= "2.15.1")  utils::globalVariables(c("vcov", "grp"))

#' Intraclass Correlation Coefficient (ICC) from Mixed Models
#'
#' This function estimates the ICC from mixed effects models
#' estimated using \pkg{lme4}.
#'
#' @param dv A character string giving the variable name of
#'   the dependent variable.
#' @param id A character vector of length one or more giving
#'   the ID variable(s).  Can be more than one.
#' @param data A data.table containing the variables
#'   used in the formula.  This is a required argument.
#'   If a data.frame, it will silently coerce to a data.table.
#'   If not a data.table or data.frame, it will attempt to coerce,
#'   with a message.
#' @param family A character vector giving the family to use
#'   for the model.  Currently only supports
#'   \dQuote{gaussian} or \dQuote{binomial}.
#' @return A data table of the ICCs
#' @references For details, see
#' Campbell, M. K., Mollison, J., and Grimshaw, J. M. (2001)
#' <doi:10.1002/1097-0258(20010215)20:3%3C391::AID-SIM800%3E3.0.CO;2-Z>
#' "Cluster trials in implementation research: estimation of intracluster correlation coefficients and sample size."
#' @keywords multivariate
#' @importFrom lme4 lmer glmer
#' @importFrom nlme VarCorr
#' @importFrom stats binomial as.formula
#' @importFrom data.table as.data.table data.table := copy
#' @export
#' @examples
#' iccMixed("mpg", "cyl", mtcars)
#' iccMixed("mpg", "cyl", data.table::as.data.table(mtcars))
#' iccMixed("mpg", "cyl", data.table::as.data.table(mtcars), family = "gaussian")
#' iccMixed("mpg", c("cyl", "am"), data.table::as.data.table(mtcars))
#' iccMixed("am", "cyl", data.table::as.data.table(mtcars), family = "binomial")
iccMixed <- function(dv, id, data, family = c("gaussian", "binomial")) {
  if (!is.data.table(data)) {
    if (is.data.frame(data)) {
      data <- as.data.table(data)
    } else {
      message("Attempting to coerce data to a data.table")
      data <- as.data.table(data)
    }
  }
  stopifnot(all(c(dv, id) %in% names(data)))
  stopifnot(is.character(dv))
  stopifnot(all(is.character(id)))
  stopifnot(identical(length(dv), 1L))
  stopifnot(length(id) >= 1L)

  d <- copy(data[, c(dv, id), with = FALSE])

  f <- sprintf("%s ~ 1 + %s", dv, paste(paste0("(1 | ", id, ")"), collapse = " + "))

  family <- match.arg(family)

  ## constant estimate of residual variance for logistic model
  ## on the 'latent variable' scale
  res.binom <- (pi^2) / 3

  m <- switch(family,
              gaussian = lmer(formula = as.formula(f), data = d, REML = TRUE),
              binomial = glmer(formula = as.formula(f), data = d, family = binomial())
              )

  est <- as.data.table(as.data.frame(VarCorr(m)))[, .(grp, vcov)]

  if (identical(family, "binomial")) {
    est <- rbind(est, est[1])
    est[nrow(est), c("grp", "vcov") := .("Residual", res.binom)]
  }

  est[, .(Var = grp, Sigma = vcov, ICC = vcov / sum(vcov))]
}


#' Estimate the effective sample size from longitudinal data
#'
#' This function estimates the (approximate) effective sample
#' size.
#'
#' @param n The number of unique/independent units of observation
#' @param k The (average) number of observations per unit
#' @param icc The estimated ICC.  If missing, will
#'   estimate (and requires that the family argument be
#'   correctly specified).
#' @param dv A character string giving the variable name of
#'   the dependent variable.
#' @param id A character vector of length one giving
#'   the ID variable.
#' @param data A data.table containing the variables
#'   used in the formula.  This is a required argument.
#'   If a data.frame, it will silently coerce to a data.table.
#'   If not a data.table or data.frame, it will attempt to coerce,
#'   with a message.
#' @param family A character vector giving the family to use
#'   for the model.  Currently only supports
#'   \dQuote{gaussian} or \dQuote{binomial}.
#' @return A data.table including the effective sample size.
#' @references For details, see
#' Campbell, M. K., Mollison, J., and Grimshaw, J. M. (2001)
#' <doi:10.1002/1097-0258(20010215)20:3%3C391::AID-SIM800%3E3.0.CO;2-Z>
#' "Cluster trials in implementation research: estimation of intracluster correlation coefficients and sample size."
#' @keywords multivariate
#' @export
#' @examples
#' ## example where n, k, and icc are estimated from the data
#' ## provided, partly using iccMixed function
#' nEffective(dv = "mpg", id = "cyl", data = mtcars)
#'
#' ## example where n, k, and icc are known (or being 'set')
#' ## useful for sensitivity analyses
#' nEffective(n = 60, k = 10, icc = .6)
nEffective <- function(n, k, icc, dv, id, data, family = c("gaussian", "binomial")) {
  if (any(missing(n), missing(k), missing(icc))) {
    if (!is.data.table(data)) {
      if (is.data.frame(data)) {
        data <- as.data.table(data)
      } else {
        message("Attempting to coerce data to a data.table")
        data <- as.data.table(data)
      }
    }
    stopifnot(all(c(dv, id) %in% names(data)))
    stopifnot(is.character(dv))
    stopifnot(all(is.character(id)))
    stopifnot(identical(length(dv), 1L))
    stopifnot(identical(length(id), 1L))

    d <- copy(data[, c(dv, id), with = FALSE])

    if (missing(icc)) {
      icc <- iccMixed(dv = dv, id = id, data = data, family = family)$ICC[1]
    }

    if (missing(n)) {
      n <- length(unique(data[[id]]))
    }

    if (missing(k)) {
      k <- nrow(data) / n
    }
  }

  neff <- (n * k) / ((1 + (k - 1) * icc))

  data.table(
    Type = c("Effective Sample Size", "Independent Units", "Total Observations"),
    N = c(neff, n, n * k))
}

#' Function to calculate the mean and deviations from mean
#'
#' Tiny helper function to calculate the mean and
#' deviations from the mean, both returned as a list.
#' Works nicely with data.table to calculate a between and
#' within variable.
#'
#' @param x A vector, appropriate for the \code{mean}
#'   function.
#' @param na.rm A logical, whether to remove missing
#'   or not.  Defaults to \code{TRUE}.
#' @return A list of the mean (first element) and deviations
#'   from the mean (second element).
#' @export
#' @examples
#' ## simple example showing what it does
#' meanDeviations(1:10)
#'
#' ## example use case, applied to a data.table
#' library(data.table)
#' d <- as.data.table(iris)
#' d[, c("BSepal.Length", "WSepal.Length") := meanDeviations(Sepal.Length),
#'   by = Species]
#' str(d)
#'
#' rm(d)
meanDeviations <- function(x, na.rm = TRUE) {
  m <- mean(x, na.rm = na.rm)
  list(m, x - m)
}

#' Estimate the autocorrelation by unit (ID)
#'
#' This function estimates the autocorrelation over time in a time
#' series by a higher level unit, given by ID.
#'
#' @param xvar A character string giving the variable name of
#'   the variable to calculate autocorrelations on.
#' @param timevar A character string giving the variable name of
#'   the time variable.
#' @param idvar A character string giving the variable name of
#'   the ID variable.  Can be missing if only one time series
#'   provided, in which case one will be created.
#' @param data A data.table containing the variables
#'   used in the formula.  This is a required argument.
#'   If a data.frame, it will silently coerce to a data.table.
#'   If not a data.table or data.frame, it will attempt to coerce,
#'   with a message.
#' @param lag.max An integer of the maximum lag to estimate. Must be
#'   equal to or greater than the number of observations
#'   for all IDs in the dataset.
#' @param na.function A character string giving the name of the function
#'   to use to address any missing data.  Functions come from the
#'   \pkg{zoo} package, and must be one of:
#'   \dQuote{na.approx}, \dQuote{na.spline}, \dQuote{na.locf}.
#' @param ... Additional arguments passed to \code{zoo}.
#' @return A data.table of the estimated autocorrelations by ID and lag
#' @keywords multivariate
#' @importFrom data.table copy is.data.table as.data.table data.table
#' @importFrom zoo zoo na.approx na.spline na.locf
#' @importFrom stats acf
#' @export
#' @examples
#' ## example 1
#' dat <- data.table::data.table(
#'   x = sin(1:30),
#'   time = 1:30,
#'   id = 1)
#' acfByID("x", "time", "id", data = dat)
#'
#' ## example 2
#' dat2 <- data.table::data.table(
#'   x = c(sin(1:30), sin((1:30)/10)),
#'   time = c(1:30, 1:30),
#'   id = rep(1:2, each = 30))
#' dat2$x[4] <- NA
#'
#' res <- acfByID("x", "time", "id", data = dat2, na.function = "na.approx")
#'
#' ggplot2::ggplot(res, ggplot2::aes(factor(Lag), AutoCorrelation)) +
#'   ggplot2::geom_boxplot()
#'
#' ## clean up
#' rm(dat, dat2, res)
acfByID <- function(xvar, timevar, idvar, data, lag.max = 10L,
                    na.function = c("na.approx", "na.spline", "na.locf"), ...) {
  if (!is.data.table(data)) {
    if (is.data.frame(data)) {
      data <- as.data.table(data)
    } else {
      message("Attempting to coerce data to a data.table")
      data <- as.data.table(data)
    }
  }

  stopifnot(is.integer(lag.max))
  stopifnot(is.character(xvar))
  stopifnot(is.character(timevar))

  stopifnot(all(c(xvar, timevar) %in% names(data)))
  stopifnot(identical(length(xvar), 1L))
  stopifnot(identical(length(timevar), 1L))

  na.function <- match.arg(na.function)
  na.function <- switch(na.function,
                        na.approx = na.approx,
                        na.spline = na.spline,
                        na.locf = na.locf)

  if (!missing(idvar)) {
    stopifnot(is.character(idvar))
    stopifnot(idvar %in% names(data))
    stopifnot(identical(length(idvar), 1L))

    d <- copy(data[, c(xvar, timevar, idvar), with = FALSE])
  } else {
    d <- copy(data[, c(xvar, timevar), with = FALSE])
    idvar <- "ID"
    while(idvar %in% names(d)) {
      idvar <- paste0("TMP_", idvar)
    }
    d[, (idvar) := 1L]
  }

  d[, .(
    Variable = xvar,
    Lag = 0:lag.max,
    AutoCorrelation = acf(na.function(zoo(get(xvar), order.by = get(timevar))),
                          lag.max = lag.max, plot = FALSE, ...)$acf[, 1, 1]),
    by = idvar]
}

#' Weighted Simple Moving Average
#'
#' This function estimates the simple moving average for a specific window
#' and weights it with a variety of optional decays (e.g., exponential, linear, none).
#' Whether to omit missing data or not is based on the missing threshold, which is a
#' proportion and indicates the tolerance. If the weighted proportion missing exceeds
#' this threshold, then that observation is missing, otherwise, missing data are excluded
#' and the weighted simple moving average calculated on the non missing data.
#'
#' @param x Time series data on which to calculate the weighted simple moving average.
#' It is assumed that these data are in the correct order and that time is
#' equally spaced. Any missing data should be filled in with NAs.
#' @param window An integer indicating the size of the window to use.
#' This window will include the current value.
#' @param decay A character string indicating the type of decay to use on the weights.
#' @param alpha An optional value. Not needed for \code{decay} = \dQuote{none}, but it
#' is required for the exponential and linear decay. For exponential and linear decay,
#' alpha should range between 0 and 1. 0 will result in no decay.
#' @param missThreshold A numeric value indicating the proportion of data that can be
#' missing for a given window before the resulting simple moving average is set to
#' missing. This is a proportion of the weighted data, so not all data points will
#' necessarily be equally weighted.
#' @return A numeric vector of the weighted simple moving averages
#' @keywords descriptives
#' @export
#' @examples
#' dweights <- expand.grid(
#'   time = 0:10,
#'   alpha = seq(0, 1, by = .1))
#'
#' library(ggplot2)
#'
#' ggplot(dweights, aes(time, (1 - alpha)^time, colour = factor(alpha))) +
#'   geom_line() + geom_point() + theme_bw() +
#'   scale_x_reverse() +
#'   theme(legend.position = "bottom") +
#'   ggtitle("Exponential Decay in Weights")
#'
#' ggplot(dweights, aes(time, pmax(1 - alpha * time, 0), colour = factor(alpha))) +
#'   geom_line() + geom_point() + theme_bw() +
#'   scale_x_reverse() +
#'   theme(legend.position = "bottom") +
#'   ggtitle("Linear Decay in Weights")
#'
#' weighted.sma(c(1, 2, 3, 4, 5),
#'              window = 3L, decay = "none",
#'              missThreshold = 0)
#'
#' weighted.sma(c(1, 2, 3, 4, 5),
#'              window = 3L, decay = "exponential",
#'              alpha = 0, missThreshold = 0)
#'
#' weighted.sma(c(1, 2, 3, 4, 5),
#'              window = 3L, decay = "linear",
#'              alpha = 0, missThreshold = 0)
#'
#' weighted.sma(c(1, 2, 3, 4, 5),
#'              window = 3L, decay = "exponential",
#'              alpha = .1, missThreshold = 0)
#'
#' weighted.sma(c(1, 2, 3, 4, 5),
#'              window = 3L, decay = "exponential",
#'              alpha = .5, missThreshold = 0)
#'
#' weighted.sma(c(1, 2, 3, 4, 5),
#'              window = 3L, decay = "linear",
#'              alpha = .1, missThreshold = 0)
#'
#' weighted.sma(c(1, 2, 3, 4, 5),
#'              window = 3L, decay = "linear",
#'              alpha = .3, missThreshold = 0)
#'
#' weighted.sma(c(1, NA, NA, 4, 5),
#'              window = 4L, decay = "exponential",
#'              alpha = .4, missThreshold = .4)
#'
#' ## clean up
#' rm(dweights)
weighted.sma <- function(x, window, decay = c("exponential", "linear", "none"), alpha, missThreshold = 0) {
  stopifnot(identical(length(window), 1L))
  if (!is.integer(window)) {
    stopifnot(as.integer(window) == window)
  }
  stopifnot(window > 0)

  n <- length(x)

  if (n < window) {
    out <- rep(NA_real_, n)
  } else {
    decay <- match.arg(decay)

    ## window minus 1 used for selecting
    window1 <- window - 1

    w <- switch(decay,
                exponential = {stopifnot(alpha >= 0); (1 - alpha)^(window1:0)},
                linear = {stopifnot(alpha >= 0); pmax(1 - alpha * (window1:0), 0)},
                none = rep(1, window))

    if (any(w < 1e-5)) {
      warning("At least one weight is effectively 0 leading to an effective window narrower than specified\nConsider decreasing alpha.")
    }

    ## adjust missing threshold proportion
    ## to absolute values based on sum of weights
    missThresholdUse <- missThreshold * sum(w)

    ## initialize logical vector of whether x is missing
    xna <- is.na(x)

    ## initialize output vector
    out <- numeric(n)
    ## set first few to missing
    out[1:window1] <- NA_real_

    for (i in window:n) {
      usex <- x[(i - window1):i]
      m <- xna[(i - window1):i]
      msum <- sum(m * w)
      if (msum > missThresholdUse) {
        out[i] <- NA_real_
      } else {
        usex <- usex[!m]
        usew <- w[!m]
        out[i] <- sum(usex * usew) / sum(usew)
      }
    }
  }

  return(out)
}

## clear R CMD CHECK notes
if (getRversion() >= "2.15.1")  utils::globalVariables(c("var", "Estimate", "LL", "UL"))

#' Summarize a Variable in a Long Dataset by ID
#' 
#' @param data A data.table object, long format
#' @param var A character string, the name of the variable to summarize
#' @param CI A numeric value, the confidence interval to use. Default is .95.
#' @param robust A logical. Default is \code{FALSE}.
#'   If \code{TRUE}, the function will use the median as the estimate.
#'   If \code{FALSE}, the function will use the mean as the estimate.
#' @param idvar A character string, the name of the grouping variable
#' @return A data.table object with the mean/median, lower limit, and upper limit of the variable
#'   specified in \code{var} for each level of the grouping variable specified in \code{idvar}.
#' @importFrom stats quantile median
#' @keywords internal
.summary.ID <- function(data, var, idvar, CI = .95, robust = FALSE) {
  stopifnot(is.data.table(data))
  stopifnot(nrow(data) > 0L)

  stopifnot(identical(length(var), 1L))
  stopifnot(identical(length(idvar), 1L))
  stopifnot(var %in% names(data))
  stopifnot(idvar %in% names(data))
  stopifnot(var != idvar)

  stopifnot(identical(length(CI), 1L))
  stopifnot(is.numeric(CI))
  stopifnot(CI > 0 & CI < 1)

  lowerlimit <- (1 - CI) / 2
  upperlimit <- 1 - lowerlimit

  stopifnot(identical(length(robust), 1L))

  if (isTRUE(robust)) {
    centerfun <- median
  } else {
    centerfun <- mean
  }

  out <- data[, .(
    Estimate = centerfun(get(var)),
    LL = quantile(get(var), lowerlimit),
    UL = quantile(get(var), upperlimit)),
    by = idvar][order(Estimate)]
  out[, (idvar) := factor(get(idvar), levels = get(idvar))]
  return(out)
}