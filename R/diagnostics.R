## clear R CMD CHECK notes
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
           "Residuals", "Index", ".N",
           "originalindex"))
}

#' residualDiagnostics methods for merMod objects
#'
#' @param object An object with class \code{merMod}. Currently only
#'   \code{lmer()} models are supported.
#' @param ev.perc The extreme value percentile to use. Defaults to
#'   .001.
#' @param robust A logical value, whether to use robust estimates or not.
#'   Defaults to \code{FALSE}.
#' @param distr A character string specifying the assumed distribution.
#'   Currently \dQuote{normal}, but may expand in the future if
#'   \code{glmer()} models are supported.
#' @param standardized A logical value whether to use standardized
#'   residual values or not. Defaults to \code{TRUE}.
#' @param ... Additional arguments. Not currently used.
#' @return A logical (\code{is.residualDiagnostics}) or
#'   a residualDiagnostics object (list) for
#'   \code{as.residualDiagnostics} and \code{residualDiagnostics}.
#' @importFrom stats model.frame resid fitted coef predict residuals
#' @importFrom data.table data.table := as.data.table
#' @importFrom JWileymisc residualDiagnostics testDistribution is.residualDiagnostics as.residualDiagnostics .quantilePercentiles
#' @method residualDiagnostics merMod
#' @export
#' @examples
#' \dontrun{
#'   sleep[1,1] <- NA
#'   m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep)
#'
#'  residualDiagnostics(m)$Residuals
#'
#' #  gm1 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
#' #    data = lme4::cbpp, family = binomial)
#' # residualDiagnostics(gm1) ## currently an error
#'
#' rm(m, sleep)
#'
#' }
residualDiagnostics.merMod <- function(object, ev.perc = .001,
                                   robust = FALSE, distr = "normal",
                                   standardized = TRUE, ...) {
  if (isGLMM(object)) stop("currently glmer() models are not supported")
  d.frame <- model.frame(object)
  naaction <- attr(d.frame, "na.action")
  if (isFALSE(is.null(naaction))) {
    if (isTRUE(inherits(naaction, "omit"))) {
      origindex <- index <- 1:(nrow(d.frame) + length(naaction))
      index[naaction] <- NA
      index[-naaction] <- 1:nrow(d.frame)
      key <- data.table(
        originalindex = origindex,
        index = index)[!is.na(index)]
    }
  } else {
    key <- data.table(
      originalindex = 1:nrow(d.frame),
      index = 1:nrow(d.frame))[!is.na(index)]
  }

  d.frame <- as.data.table(d.frame)
  dv <- names(d.frame)[1]

  d.res <- data.table(
    Residuals = residuals(object, type = "response", scaled = standardized),
    Predicted = fitted(object))

  d.hat <- .quantilePercentiles(
    data = d.res,
    LL = .1, UL = .9)

  d.dist <- testDistribution(
    x = d.res$Residuals,
    distr = "normal",
    na.rm = TRUE,
    extremevalues = "theoretical",
    ev.perc = ev.perc,
    use = "complete.obs",
    robust = robust)

  d.res[!is.na(Residuals), isEV := d.dist$Data[order(OriginalOrder), isEV]]
  d.res[, Index := 1:.N]

  ## fix the index to match original data if missing data existed and were omitted
  if (isFALSE(is.null(naaction))) {
    if (isTRUE(inherits(naaction, "omit"))) {
      d.res[, Index := key[, originalindex]]
      d.dist$Data[, OriginalOrder := key[match(OriginalOrder, index), originalindex]]
    }
  }

  out <- list(na.omit(d.res), d.frame, d.hat, d.dist, dv)
  attr(out, "augmentClass") <- "merMod"

  as.residualDiagnostics(out)
}


## clear R CMD CHECK notes
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".SD", "isEV", "OriginalOrder"))

#' modelDiagnostics method for merMod objects
#'
#' This function creates a number of diagnostics
#' for \code{merMod} models.
#'
#' @param object A fitted model object, either of class merMod from
#'   the lme4 package or merModLmerTest from the lmerTest package.
#' @param ev.perc A real number between 0 and 1 indicating the
#'   proportion of the theoretical distribution beyond which
#'   values are considered extreme values (possible outliers).
#'   Defaults to .001.
#' @param robust Whether to use robust mean and standard deviation estimates
#'   for normal distribution
#' @param distr A character string given the assumed distribution.
#'   Passed on to \code{\link{testDistribution}}.
#'   Defaults to \dQuote{normal}.
#' @param standardized A logical whether to use standardized residuals.
#'   Defaults to \code{TRUE} generally where possible but may depend on
#'   method.
#' @param ... Additional arguments, not currently used.
#' @return A logical (\code{is.modelDiagnostics}) or
#'   a modelDiagnostics object (list) for
#'   \code{as.modelDiagnostics} and \code{modelDiagnostics}.
#' @importFrom JWileymisc as.na testDistribution modelDiagnostics as.modelDiagnostics
#' @importFrom nlme VarCorr ranef
#' @importFrom data.table as.data.table :=
#' @importFrom stats na.omit
#' @method modelDiagnostics merMod
#' @keywords plot
#' @export
#' @examples
#' \dontrun{
#'   sleep[1,1] <- NA
#'   m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep)
#'
#'  md <- modelDiagnostics(m, ev.perc = .1)
#'  md$extremeValues
#' class(md)
#'
#' plot(md)
#'
#' data(aces_daily, package = "JWileymisc")
#' m <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID), data = aces_daily)
#' md <- modelDiagnostics(m, ev.perc = .1)
#'
#' #  gm1 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
#' #    data = lme4::cbpp, family = binomial)
#' # modelDiagnostics(gm1) ## currently an error
#'
#' rm(m, md, sleep)
#'
#' }
modelDiagnostics.merMod <- function(object, ev.perc = .001,
                                   robust = FALSE, distr = "normal",
                                   standardized = TRUE, ...) {
  if (isGLMM(object)) stop("currently glmer() models are not supported")

  x <- residualDiagnostics(object,
                           ev.perc = ev.perc,
                           robust = robust,
                           distr = distr,
                           standardized = standardized)

  idvars <- names(VarCorr(object))

  ## data for outliers
  d.extreme <- x$Frame[1, c(x$Outcome, idvars),
                       with = FALSE][, lapply(.SD, as.na)]
  d.extreme[, OriginalOrder := NA_integer_]

  if ("EffectType" %in% names(d.extreme)) {
    stop("EffectType is used internally and cannot be a variable or ID in the model")
  }
  d.extreme[, EffectType := NA_character_]

  x$Frame[, isEV := x$testDistribution$Data[order(OriginalOrder), isEV]]
  x$Frame[, OriginalOrder := x$testDistribution$Data[order(OriginalOrder), OriginalOrder]]

  if (isTRUE(any(x$Frame[, isEV] == "Yes"))) {
    d.extreme <- rbind(d.extreme,
                       cbind(x$Frame[isEV == "Yes", c(x$Outcome, idvars, "OriginalOrder"),
                                     with = FALSE],
                         EffectType = "Residuals"))
  }

  p.ranef <- list()
  for (n in idvars) {
    tmp <- subset(coef(object)[[n]],
                  select = names(ranef(object)[[n]]))
    for (n2 in names(tmp)) {
      p.tmpranef <- testDistribution(tmp[[n2]],
                              varlab = "Random Effects",
                              plot = FALSE, extremevalues = "theoretical",
                              ev.perc = ev.perc)
      p.ranef <- c(p.ranef,
                   list(list(
                     Title = paste(n, ":", n2),
                     p.tmpranef)))
      if (any(p.tmpranef$Data[, isEV] == "Yes")) {
        d.extreme <- rbind(d.extreme,
                           cbind(x$Frame[x$Frame[[n]] %in% rownames(tmp)[
                                                             p.tmpranef$Data[isEV == "Yes", OriginalOrder]],
                                         c(x$Outcome, idvars, "OriginalOrder"), with = FALSE],
                                 EffectType = paste("Random Effect", n, ":", n2)))
      }
    }
    if (ncol(tmp) > 1) {
      p.tmpranef <- testDistribution(tmp, distr = "mvnormal",
                              varlab = "Random Effects",
                              plot = FALSE, extremevalues = "theoretical",
                              ev.perc = ev.perc)
      p.ranef <- c(p.ranef,
                   list(list(
                     Title = paste(n, ":", "MV Normal"),
                     p.tmpranef)))
      if (any(p.tmpranef$Data[, isEV] == "Yes")) {
        d.extreme <- rbind(d.extreme,
                           cbind(x$Frame[x$Frame[[n]] %in% rownames(tmp)[
                                                             p.tmpranef$Data[isEV == "Yes", OriginalOrder]],
                                         c(x$Outcome, idvars, "OriginalOrder"), with = FALSE],
                                 EffectType = paste("Multivariate Random Effect", n)))
      }
    }
  }

  setnames(d.extreme, old = "OriginalOrder", new = "Index")

  out <- list(x, p.ranef, na.omit(d.extreme))
  attr(out, "augmentClass") <- "merMod"

  as.modelDiagnostics(out)
}




#' residualDiagnostics methods for lme objects
#'
#' @param object An object with class \code{lme}.
#' @param ev.perc The extreme value percentile to use. Defaults to
#'   .001.
#' @param robust A logical value, whether to use robust estimates or not.
#'   Defaults to \code{FALSE}.
#' @param distr A character string specifying the assumed distribution.
#'   Currently \dQuote{normal}, but future options may be supported in the
#'   future.
#' @param standardized A logical value whether to use standardized pearson
#'   residual values or not. Defaults to \code{TRUE}.
#' @param ... Additional arguments. Not currently used.
#' @return A logical (\code{is.residualDiagnostics}) or
#'   a residualDiagnostics object (list) for
#'   \code{as.residualDiagnostics} and \code{residualDiagnostics}.
#' @importFrom nlme getGroupsFormula
#' @importFrom stats get_all_vars
#' @importFrom data.table data.table := as.data.table
#' @importFrom JWileymisc residualDiagnostics testDistribution is.residualDiagnostics as.residualDiagnostics .quantilePercentiles
#' @method residualDiagnostics lme
#' @export
#' @examples
#' \dontrun{
#'   sleep[1,1] <- NA
#'   m <- nlme::lme(extra ~ group, data = sleep, random = ~ 1 | ID,
#'     na.action = na.omit)
#'
#'  residualDiagnostics(m)$Residuals
#'
#'   m <- nlme::lme(hp ~ mpg, data = mtcars, random = ~ 1 | cyl,
#'     na.action = na.omit)
#'  residualDiagnostics(m)$Residuals
#'
#' rm(m, sleep)
#'
#' }
residualDiagnostics.lme <- function(object, ev.perc = .001,
                                   robust = FALSE, distr = "normal",
                                   standardized = TRUE, ...) {
  d.frame <- na.omit(cbind(
    get_all_vars(formula(object), data = object$data),
    get_all_vars(getGroupsFormula(object), data = object$data)))

  naaction <- object$na.action
  if (isFALSE(is.null(naaction))) {
    if (isTRUE(inherits(naaction, "omit"))) {
      origindex <- index <- 1:(nrow(d.frame) + length(naaction))
      index[naaction] <- NA
      index[-naaction] <- 1:nrow(d.frame)
      key <- data.table(
        originalindex = origindex,
        index = index)[!is.na(index)]
    }
  } else {
    key <- data.table(
      originalindex = 1:nrow(d.frame),
      index = 1:nrow(d.frame))[!is.na(index)]
  }

  d.frame <- as.data.table(d.frame)
  dv <- names(d.frame)[1]

  d.res <- data.table(
    Residuals = residuals(object, type = if (standardized) "pearson" else "response"),
    Predicted = fitted(object))

  d.hat <- .quantilePercentiles(
    data = d.res,
    LL = .1, UL = .9)

  d.dist <- testDistribution(
    x = d.res$Residuals,
    distr = "normal",
    na.rm = TRUE,
    extremevalues = "theoretical",
    ev.perc = ev.perc,
    use = "complete.obs",
    robust = robust)

  d.res[!is.na(Residuals), isEV := d.dist$Data[order(OriginalOrder), isEV]]
  d.res[, Index := 1:.N]

  ## fix the index to match original data if missing data existed and were omitted
  if (isFALSE(is.null(naaction))) {
    if (isTRUE(inherits(naaction, "omit"))) {
      d.res[, Index := key[, originalindex]]
      d.dist$Data[, OriginalOrder := key[match(OriginalOrder, index), originalindex]]
    }
  }

  out <- list(na.omit(d.res), d.frame, d.hat, d.dist, dv)
  attr(out, "augmentClass") <- "lme"

  as.residualDiagnostics(out)
}

#' modelDiagnostics method for lme objects
#'
#' This function creates a number of diagnostics
#' for \code{lme} models.
#'
#' @param object A fitted model object  of class lme.
#' @param ev.perc A real number between 0 and 1 indicating the
#'   proportion of the theoretical distribution beyond which
#'   values are considered extreme values (possible outliers).
#'   Defaults to .001.
#' @param robust Whether to use robust mean and standard deviation estimates
#'   for normal distribution
#' @param distr A character string given the assumed distribution.
#'   Passed on to \code{\link{testDistribution}}.
#'   Defaults to \dQuote{normal}.
#' @param standardized A logical whether to use standardized pearson residuals.
#'   Defaults to \code{TRUE} generally where possible but may depend on
#'   method.
#' @param ... Additional arguments, not currently used.
#' @return A logical (\code{is.modelDiagnostics}) or
#'   a modelDiagnostics object (list) for
#'   \code{as.modelDiagnostics} and \code{modelDiagnostics}.
#' @importFrom JWileymisc as.na testDistribution modelDiagnostics as.modelDiagnostics
#' @importFrom nlme VarCorr ranef getGroupsFormula random.effects
#' @importFrom data.table as.data.table :=
#' @importFrom stats na.omit get_all_vars
#' @method modelDiagnostics lme
#' @keywords plot
#' @export
#' @examples
#' \dontrun{
#'   sleep[1,1] <- NA
#'   m <- nlme::lme(extra ~ group, data = sleep,
#'    random = ~ 1 | ID, na.action = "na.omit")
#'
#'  md <- modelDiagnostics(m, ev.perc = .1)
#'  md$extremeValues
#'
#'  plot(md)
#'
#' data(aces_daily, package = "JWileymisc")
#' m <- nlme::lme(PosAff ~ STRESS, data = aces_daily,
#'   random = ~ 1 + STRESS | UserID, na.action = "na.omit")
#' md <- modelDiagnostics(m, ev.perc = .001)
#' md$extremeValues
#' plot(md)
#'
#'
#' m <- nlme::lme(extra ~ 1, data = sleep, random = ~ 1 | ID/group,
#'   na.action = "na.omit")
#'
#' md <- modelDiagnostics(m, ev.perc = .1)
#' md$extremeValues
#' plot(md)

#' rm(m, md, sleep)
#'
#' }
modelDiagnostics.lme <- function(object, ev.perc = .001,
                                   robust = FALSE, distr = "normal",
                                   standardized = TRUE, ...) {

  x <- residualDiagnostics(object,
                           ev.perc = ev.perc,
                           robust = robust,
                           distr = distr,
                           standardized = standardized)

  idvars <- names(get_all_vars(nlme::getGroupsFormula(object), data = object$data))

  ## data for outliers
  d.extreme <- x$Frame[1, c(x$Outcome, idvars),
                       with = FALSE][, lapply(.SD, as.na)]
  d.extreme[, OriginalOrder := NA_integer_]

  if ("EffectType" %in% names(d.extreme)) {
    stop("EffectType is used internally and cannot be a variable or ID in the model")
  }
  d.extreme[, EffectType := NA_character_]

  x$Frame[, isEV := x$testDistribution$Data[order(OriginalOrder), isEV]]
  x$Frame[, OriginalOrder := x$testDistribution$Data[order(OriginalOrder), OriginalOrder]]

  if (isTRUE(any(x$Frame[, isEV] == "Yes"))) {
    d.extreme <- rbind(d.extreme,
                       cbind(x$Frame[isEV == "Yes", c(x$Outcome, idvars, "OriginalOrder"),
                                     with = FALSE],
                         EffectType = "Residuals"))
  }

  p.ranef <- list()
  for (n in idvars) {
    if (identical(length(idvars), 1L)) {
      tmp <- random.effects(object)
    } else {
      tmp <- random.effects(object)[[n]]
    }
    fe <- fixef(object)
    for (n2 in names(tmp)) {
      if (isTRUE(n2 %in% names(fe))) {
        tmp[[n2]] <- tmp[[n2]] + fe[[n2]]
      }
      p.tmpranef <- testDistribution(tmp[[n2]],
                              varlab = "Random Effects",
                              plot = FALSE, extremevalues = "theoretical",
                              ev.perc = ev.perc)
      p.ranef <- c(p.ranef,
                   list(list(
                     Title = paste(n, ":", n2),
                     p.tmpranef)))
      if (any(p.tmpranef$Data[, isEV] == "Yes")) {
        d.extreme <- rbind(d.extreme,
                           cbind(x$Frame[x$Frame[[n]] %in% rownames(tmp)[
                                                             p.tmpranef$Data[isEV == "Yes", OriginalOrder]],
                                         c(x$Outcome, idvars, "OriginalOrder"), with = FALSE],
                                 EffectType = paste("Random Effect", n, ":", n2)))
      }
    }
    if (ncol(tmp) > 1) {
      p.tmpranef <- testDistribution(tmp, distr = "mvnormal",
                              varlab = "Random Effects",
                              plot = FALSE, extremevalues = "theoretical",
                              ev.perc = ev.perc)
      p.ranef <- c(p.ranef,
                   list(list(
                     Title = paste(n, ":", "MV Normal"),
                     p.tmpranef)))
      if (any(p.tmpranef$Data[, isEV] == "Yes")) {
        d.extreme <- rbind(d.extreme,
                           cbind(x$Frame[x$Frame[[n]] %in% rownames(tmp)[
                                                             p.tmpranef$Data[isEV == "Yes", OriginalOrder]],
                                         c(x$Outcome, idvars, "OriginalOrder"), with = FALSE],
                                 EffectType = paste("Multivariate Random Effect", n)))
      }
    }
  }

  setnames(d.extreme, old = "OriginalOrder", new = "Index")

  out <- list(x, p.ranef, na.omit(d.extreme))
  attr(out, "augmentClass") <- "lme"

  as.modelDiagnostics(out)
}


