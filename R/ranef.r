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


## clear R CMD CHECK notes
if (getRversion() >= "2.15.1")  utils::globalVariables(c("iteration"))

#' @rdname ranef2long
#' @param d A \code{ranef} object
#' @param i an integer, which random effect to pull out
#' @importFrom data.table as.data.table melt
.re.data <- function(d, i, idvar) {
  xw <- as.data.table(t(d[, , i]))
  xw[, (idvar) := dimnames(d)[[2]]]
  xlong <- melt(xw, id.vars = idvar,
    value.name = dimnames(d)[[3]][i],
    variable.name = "iteration")
  xlong[, iteration := as.integer(iteration)]
  return(xlong)
}

## clear R CMD CHECK notes
if (getRversion() >= "2.15.1")  utils::globalVariables(c("e", ".x", "A", "B",
  "interceptA", "interceptB", "sigmaA", "sigmaB", "x", "y", "vars"))

#' Create data and plots for \code{\link[brms]{brm}} random effect models
#' 
#' @param object a \code{\link[brms]{brmsfit-class}} object
#' @param usevars a character vector of random effects to plot
#' @param newdata a \code{data.table} object with the data used to generate the random effects, this is used as an anchor for the random intercepts so they have a meaningful 0 point
#' @param idvar a character string specifying the grouping variable name for the random effects 
#' @param CI a numeric value between 0 and 1 specifying the interval to use. Defaults to 0.95.
#' @param robust a logical value indicating whether to use robust estimates or not. Defaults to FALSE.
#'   Passed on to \code{\link[brms]{posterior_summary}} and \code{\link{.summary.ID}}.
#' @return a list with the following components:
#' * \code{plot}: a list of ggplot objects
#' * \code{plotdat}: a list of data.table objects with the data used to generate the plots
#' * \code{relong}: a data.table object with the random effects in long format
#' * \code{yhat}: a list of data.table objects with the expected values for the random effects
#' * \code{usevars}: a character vector of the random effects to plot
#' * \code{idvar}: a character string specifying the grouping variable name for the random effects
#' @importFrom nlme ranef fixef
#' @importFrom brms posterior_epred posterior_summary
#' @importFrom testthat expect_true
#' @importFrom extraoperators %ain%
#' @importFrom ggplot2 ggplot aes annotate geom_hex geom_pointrange geom_hline ggtitle
#' @importFrom ggplot2 stat_smooth scale_fill_continuous theme coord_flip
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous element_blank xlab ylab
#' @importFrom scales breaks_log log_trans label_math
#' @importFrom ggpubr theme_pubr
#' @importFrom data.table setkey data.table
#' @export
#' @examples
#' \donttest{
#' if (FALSE) {
#' library(data.table)
#' library(brms)
#' library(ggpubr)
#' 
#' current.seed <- .Random.seed
#' set.seed(12345)
#' nGroups <- 100
#' nObs <- 20
#' theta.location <- matrix(rnorm(nGroups * 2), nrow = nGroups, ncol = 2)
#' theta.location[, 1] <- theta.location[, 1] - mean(theta.location[, 1])
#' theta.location[, 2] <- theta.location[, 2] - mean(theta.location[, 2])
#' theta.location[, 1] <- theta.location[, 1] / sd(theta.location[, 1])
#' theta.location[, 2] <- theta.location[, 2] / sd(theta.location[, 2])
#' theta.location <- theta.location %*% chol(matrix(c(1.5, -.25, -.25, .5^2), 2))
#' theta.location[, 1] <- theta.location[, 1] - 2.5
#' theta.location[, 2] <- theta.location[, 2] + 1
#' 
#' dmixed <- data.table(
#'   x = rep(rep(0:1, each = nObs / 2), times = nGroups))
#'   dmixed[, ID := rep(seq_len(nGroups), each = nObs)]
#' 
#'   for (i in seq_len(nGroups)) {
#'     dmixed[ID == i, y := rnorm(
#'       n = nObs,
#'       mean = theta.location[i, 1] + theta.location[i, 2] * x,
#'       sd = exp(1 + theta.location[i, 1] + theta.location[i, 2] * x))
#'       ]
#'   }
#' 
#' ## note this model takes several minutes, even on a high performance machine
#' ls.me <- brm(bf(
#'   y ~ 1 + x + (1 + x | p | ID),
#'   sigma ~ 1 + x + (1 + x | p | ID)),
#'   family = "gaussian",
#'   data = dmixed, seed = 1234,
#'   silent = 2, refresh = 0, iter = 2000, warmup = 1000, thin = 1,
#'   chains = 4L, cores = 4L)
#' 
#' out <- ranefdata(
#'   ls.me,
#'   newdata = data.table(ID = dmixed$ID[1], x = 0), 
#'   usevars = c("Intercept", "x", "sigma_Intercept", "sigma_x"),
#'   idvar = "ID")
#' 
#' do.call(ggarrange, c(out$replots, ncol=2,nrow=2))
#' do.call(ggarrange, c(out$scatterplots, ncol=2,nrow=3))
#' 
#' ## set seed back to what it was
#' set.seed(current.seed)
#' 
#' ## cleanup
#' rm(current.seed, nGroups, nObs, theta.location, dmixed, ls.me, out)
#' }
#' }
ranefdata <- function(object, usevars, newdata, idvar, CI = .95, robust = FALSE) {
  if (!isTRUE(inherits(object, "brmsfit"))) {
    stop("object must be a brmsfit object")
  }

  stopifnot(identical(length(CI), 1L))
  stopifnot(is.numeric(CI))
  stopifnot(CI > 0 & CI < 1)

  lowerlimit <- (1 - CI) / 2
  upperlimit <- 1 - lowerlimit
  lowerlabel <- paste0("Q", lowerlimit * 100)
  upperlabel <- paste0("Q", upperlimit * 100)

  stopifnot(identical(length(robust), 1L))

  intercept <- grepl("Intercept", usevars)

  dpars <- family(object)$dpars
  ## only count some as "sigma" params if sigma is in dpars for the model
  ## family, this is to prevent catching variables named sigma
  if (isTRUE("sigma" %in% dpars)) {
    sigma <- grepl("^sigma.*$", usevars)
  } else {
    sigma <- rep(FALSE, length(usevars))
  }

  fes <- fixef(object, robust = robust, probs = c(lowerlimit, upperlimit))
  fes <- cbind(vars = rownames(fes), as.data.table(fes))
  setkey(fes, vars)
  setnames(fes, c(lowerlabel, upperlabel), c("LL", "UL"))

  stopifnot(usevars[!intercept] %ain% fes[, vars])

  res <- ranef(object, summary = FALSE)[[idvar]]
  relong <- ranef2long(object, idvar)

  stopifnot(usevars %ain% names(relong)) ## all usevars in the randome effects

  yhat <- vector("list", length(usevars))
  names(yhat) <- usevars

  for (i in seq_along(usevars)) {
    if (isTRUE(intercept[i])) {
      if (isTRUE(sigma[i])) {
        yhat[[i]] <- as.data.table(posterior_summary(posterior_epred(
          object,
          newdata = newdata, re_formula = NA,
          dpar = "sigma"
        ), probs = c(lowerlimit, upperlimit), robust = robust))
        setnames(yhat[[i]], c(lowerlabel, upperlabel), c("LL", "UL"))
        relong[, (usevars[i]) := exp(
          get(usevars[i]) +
            log(yhat[[i]][1, Estimate])
        )]
      } else if (isFALSE(sigma[i])) {
        yhat[[i]] <- as.data.table(posterior_summary(posterior_epred(
          object,
          newdata = newdata, re_formula = NA,
          dpar = NULL
        ), probs = c(lowerlimit, upperlimit), robust = robust))
        setnames(yhat[[i]], c(lowerlabel, upperlabel), c("LL", "UL"))
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
    plotdat[[i]] <- .summary.ID(relong, usevars[i], idvar = idvar, CI = CI, robust = robust)
    tmpplot <- ggplot(plotdat[[i]], aes(ID, Estimate, ymin = LL, ymax = UL)) +
      annotate("rect",
        xmin = -Inf, xmax = Inf,
        ymin = yhat[[i]][1, LL], ymax = yhat[[i]][1, UL],
        fill = "grey80"
      ) +
      geom_hline(yintercept = yhat[[i]][1, Estimate], linetype = "dashed") +
      geom_pointrange()

    if (sigma[i] && intercept[i]) {
      tmpplot <- tmpplot +
        scale_y_continuous(
          trans = log_trans(),
          breaks = breaks_log(n = 10, base = exp(1)),
          labels = label_math(e^.x, format = log)
        )
    }
    tmpplot <- tmpplot +
      xlab("ID") +
      ylab(usevars[i]) +
      theme_pubr() +
      theme(axis.text.y = element_blank()) +
      coord_flip()

    plot[[i]] <- tmpplot
  }

  ## make scatter plots for the random effects
  ## we only want these for all unique pairs of random effects

  ## make all pairs of usevars
  tmp <- as.data.table(expand.grid(
    A = usevars, B = usevars,
    stringsAsFactors = FALSE
  ))

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
      stat_smooth(
        method = "lm", formula = y ~ x, se = FALSE,
        colour = "white", linewidth = 2
      ) +
      scale_fill_continuous(type = "viridis") +
      theme_pubr() +
      xlab(tmp[i, A]) +
      ylab(tmp[i, B])

    if (tmp[i, interceptA] & tmp[i, sigmaA]) {
      tmpplot <- tmpplot +
        scale_x_continuous(
          trans = log_trans(),
          breaks = breaks_log(n = 10, base = exp(1)),
          labels = label_math(e^.x, format = log)
        )
    }
    if (tmp[i, interceptB] & tmp[i, sigmaB]) {
      tmpplot <- tmpplot +
        scale_y_continuous(
          trans = log_trans(),
          breaks = breaks_log(n = 10, base = exp(1)),
          labels = label_math(e^.x, format = log)
        )
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