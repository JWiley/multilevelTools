#' Convert ranef() output to long format
#' @rdname ranef2long
#' 
#' @param x A \code{brmsfit} object
#' @param idvar A character string specifying the grouping variable name for the 
#'   random effects.
#' @importFrom nlme ranef
#' @return A data.table object with the random effects in long format.
#' @keywords internal
ranef2long <- function(x, idvar) {
  x <- ranef(x, summary = FALSE)[[idvar]]
  for (i in seq_along(dimnames(x)[[3]])) {
    if (i == 1) {
      out <- .re.data(x, i, idvar)
    } else {
      tmp <- .re.data(x, i, idvar)
      out <- merge(out, tmp, by = c(idvar, "iteration"))
    }
  }
  return(out)
}

#' @rdname ranef2long
#' @param d A \code{ranef} object
#' @param i an integer, which random effect to pull out
.re.data <- function(d, i, idvar) {
  xw <- as.data.table(t(d[, , i]))
  xw[, (idvar) := dimnames(d)[[2]]]
  xlong <- melt(xw, id.vars = idvar, 
    value.name = dimnames(d)[[3]][i], 
    variable.name = "iteration")
  xlong[, iteration := as.integer(iteration)]
  return(xlong)
}

#' @name Convert re.data() output to a dataset for plotting
#' 
#' @param data A data.table object, long format random effects
#' @param var A character string, the name of the random effect to plot
#' @return A data.table object with the random effect in a format suitable for plotting
#' @keywords internal
.make.replotdat <- function(data, var) {
  out <- data[, .(
    Estimate = mean(get(var)),
    Q2.5 = quantile(get(var), .025),
    Q97.5 = quantile(get(var), .975)),
    by = ID][order(Estimate)]
  out[, ID := factor(ID, levels = ID)]
  return(out)
}

#' Create data and plots for brms random effect models
#' 
#' @param object a \code{brmsfit} objectx
#' @param usevars a character vector of random effects to plot
#' @param newdata a data.table object with the data used to generate the random effects, this is used as an anchor for the random intercepts so they have a meaningful 0 point
#' @param idvar a character string specifying the grouping variable name for the random effects 
#' @return a list with the following components:
#' \itemize{
#'   \item{plot}{a list of ggplot objects}
#'   \item{plotdat}{a list of data.table objects with the data used to generate the plots}
#'   \item{relong}{a data.table object with the random effects in long format}
#'   \item{yhat}{a list of data.table objects with the expected values for the random effects}
#'   \item{usevars}{a character vector of the random effects to plot}
#'   \item{idvar}{a character string specifying the grouping variable name for the random effects}
#' }
#' @importFrom nlme ranef fixef
#' @importFrom brms posterior_epred
#' @importFrom testthat expect_true
#' @importFrom extraoperators %ain%
#' @importFrom ggplot2 ggplot aes annotate geom_hex geom_pointrange geom_hline ggtitle
#' @importFrom ggplot2 stat_smooth scale_fill_continuous theme_classic theme axis.text.y coord_flip
#' @export
ranefdata <- function(object, usevars, newdata, idvar) {
  intercept <- grepl("Intercept", usevars)

  dpars <- family(object)$dpars
  ## only count some as "sigma" params if sigma is in dpars for the model
  ## family, this is to prevent catching variables named sigma
  if (isTRUE("sigma" %in% dpars)) {
    sigma <- grepl("^sigma.*$", usevars)
  } else {
    sigma <- rep(FALSE, length(usevars))
  }

  fes <- fixef(object)
  fes <- cbind(vars = rownames(fes), as.data.table(fes))
  setkey(fes, vars)

  expect_true(usevars[!intercept] %ain% fes[, vars])

  res <- ranef(object, summary = FALSE)[[idvar]]
  relong <- ranef2long(object, idvar)

  expect_true(usevars %ain% names(relong)) ## all usevars in the randome effects

  yhat <- vector("list", length(usevars))
  names(yhat) <- usevars

  for (i in seq_along(usevars)) {
    if (isTRUE(intercept[i])) {
      if (isTRUE(sigma[i])) {
        yhat[[i]] <- as.data.table(posterior_summary(posterior_epred(
          object,
          newdata = newdata, re_formula = NA,
          dpar = "sigma"
        )))
        relong[, (usevars[i]) := exp(
          get(usevars[i]) +
            log(yhat[[i]][1, Estimate])
        )]
      } else if (isFALSE(sigma[i])) {
        yhat[[i]] <- as.data.table(posterior_summary(posterior_epred(
          object,
          newdata = newdata, re_formula = NA,
          dpar = NULL
        )))
        relong[, (usevars[i]) := (
          get(usevars[i]) +
            (yhat[[i]][1, Estimate]))]
      }
    } else {
      yhat[[i]] <- fes[vars == usevars[i]]
      yhat[[i]][, vars := NULL]
      relong[, (usevars[i]) := (
        get(usevars[i]) +
          (yhat[[i]][1, Estimate]))]
    }
  }

  plot <- plotdat <- vector("list", length(usevars))
  names(plot) <- names(plotdat) <- usevars

  for (i in seq_along(usevars)) {
    plotdat[[i]] <- .make.replotdat(relong, usevars[i])
    tmpplot <- ggplot(plotdat[[i]], aes(ID, Estimate, ymin = Q2.5, ymax = Q97.5)) +
      annotate("rect",
        xmin = -Inf, xmax = Inf,
        ymin = yhat[[i]][1, Q2.5], ymax = yhat[[i]][1, Q97.5],
        fill = "grey80"
      ) +
      geom_hline(yintercept = yhat[[i]][1, Estimate], linetype = "dashed") +
      geom_pointrange()
      
    if (sigma[i] & intercept[i]) {
      tmpplot <- tmpplot + 
        scale_y_continuous(trans = log_trans(),
                  breaks = breaks_log(n = 10, base = exp(1)),
                  labels = label_math(e^.x, format = log))
    }
    tmpplot <- tmpplot +
      xlab("ID") +
      ylab(usevars[i]) +
      theme_classic() +
      theme(axis.text.y = element_blank()) +
      coord_flip()

    plot[[i]] <- tmpplot
  }


  ## make scatter plots for the random effects
  ## we only want these for all unique pairs of random effects

  ## make all pairs of usevars
  tmp <- as.data.table(expand.grid(
    A = usevars, B = usevars,
    stringsAsFactors = FALSE))

  ## remove any pairs that are duplicated or only reversed
  tmp <- tmp[duplicated(apply(tmp, 1, function(x) {
    paste(sort(x), collapse = "")
  }))][A != B]

  tmp[, interceptA := grepl("Intercept", A)]
  tmp[, interceptB := grepl("Intercept", B)]

  ## only count some as "sigma" params if sigma is in dpars for the model
  ## family, this is to prevent catching variables named sigma
  if (isTRUE("sigma" %in% dpars)) {
    tmp[, sigmaA := grepl("^sigma.*$", A)]
    tmp[, sigmaB := grepl("^sigma.*$", B)]
  } else {
    tmp[, sigmaA := rep(FALSE, .N)]
    tmp[, sigmaB := rep(FALSE, .N)]
  }

  scatterplot <- vector("list", nrow(tmp))

  for (i in seq_len(nrow(tmp))) {
    tmpdat <- copy(relong[, .(x = get(tmp[i, A]), y = get(tmp[i, B]))])

    tmpplot <- ggplot(tmpdat, aes(x = x, y = y)) +
      geom_hex(show.legend = FALSE) +
      stat_smooth(method = "lm", formula = y ~ x, se = FALSE,
        colour = "white", linewidth = 2) + 
      scale_fill_continuous(type = "viridis") + 
      theme_classic() + 
      xlab(tmp[i, A]) + ylab(tmp[i, B])

    if (tmp[i, interceptA] & tmp[i, sigmaA]) {
      tmpplot <- tmpplot + 
        scale_x_continuous(trans = log_trans(),
                  breaks = breaks_log(n = 10, base = exp(1)),
                  labels = label_math(e^.x, format = log))
    }
    if (tmp[i, interceptB] & tmp[i, sigmaB]) {
      tmpplot <- tmpplot + 
        scale_y_continuous(trans = log_trans(),
                  breaks = breaks_log(n = 10, base = exp(1)),
                  labels = label_math(e^.x, format = log))
    }
    scatterplot[[i]] <- tmpplot
  }

  out <- list(
    replots = plot,
    scatterplots = scatterplot,
    replotdat = plotdat,
    relong = relong,
    yhat = yhat,
    usevars = usevars,
    idvar = idvar
  )
  return(out)
}


if (FALSE) {
library(knitr)
library(data.table)
library(brms)
library(cmdstanr)
library(mice)
library(miceadds)
library(ggplot2)
library(bayesplot)
library(testthat)

dmixed <- withr::with_seed(
  seed = 12345, code = {
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
    copy(dmixed)
  })


ls.me <- brm(bf(
  y ~ 1 + x + (1 + x | p | ID),
  sigma ~ 1 + x + (1 + x | p | ID)),
  family = "gaussian",
  data = dmixed, seed = 1234,
  silent = 2, refresh = 0, iter = 4000, warmup = 1000, thin = 3,
  chains = 4L, cores = 4L, backend = "cmdstanr")

out <- ranefdata(
  ls.me,
  newdata = data.table(ID = dmixed$ID[1], x = 0), 
  usevars = c("Intercept", "x", "sigma_Intercept", "sigma_x"),
  idvar = "ID")

do.call(ggarrange, c(out$replots, ncol=2,nrow=2))
do.call(ggarrange, c(out$scatterplots, ncol=2,nrow=3))

}