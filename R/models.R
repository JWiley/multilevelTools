#' Calculate multilevel omega reliability
#'
#' This function uses multilevel structural equation modelling
#' to calculate between and within reliability using coefficient
#' omega.
#'
#' @param items A character vector giving the variables that map
#'  to the items in the scale. Note that these should be reverse
#'  scored prior to running this function.
#' @param id A character string giving the name of the variable that
#'  indicates which rows of the dataset belong to the same person
#'  or group for the multilevel analysis.
#' @param data A data table or data frame to be used for analysis.
#' @param savemodel A logical value indicating whether the underlying model
#'  should be saved and returned. Defaults to \code{FALSE}.
#' @return a list with two elements, the first, \dQuote{Results} contains the
#'  estimates for coefficient omega at the within and between level. The
#'  next element, \dQuote{Fit} contains the entire fitted model from lavaan, if
#'  \code{savemodel = TRUE}.
#' @references Geldhof, G. J., Preacher, K. J., & Zyphur, M. J. (2014)
#' <doi:10.1037/a0032138>
#' "Reliability estimation in a multilevel confirmatory factor analysis framework"
#' @export
#' @importFrom lavaan sem parameterEstimates
#' @examples
#'
#'   data(aces_daily, package = "JWileymisc")
#'   omegaSEM(
#'     items = c("COPEPrb", "COPEPrc", "COPEExp"),
#'     id = "UserID",
#'     data = aces_daily,
#'     savemodel = FALSE)
omegaSEM <- function(items, id, data, savemodel = FALSE) {
  if (length(items) < 2) {
    stop("omega requires at least two items")
  }

  llabels.within <- paste0("wl", seq_along(items))
  rlabels.within <- paste0("wr", seq_along(items))
  constraints.within <- paste(
    sprintf("%s > 0", rlabels.within),
    collapse = "\n")
  loadings.within <- paste(c(
    sprintf("NA * %s", items[[1]]),
    sprintf("%s * %s", llabels.within, items)),
    collapse = " + ")
  variances.within <- paste(sprintf(
    "%s~~%s*%s", items, rlabels.within, items),
    collapse = "\n")

  llabels.between <- paste0("bl", seq_along(items))
  rlabels.between <- paste0("br", seq_along(items))
  constraints.between <- paste(
    sprintf("%s > 0", rlabels.between),
    collapse = "\n")
  loadings.between <- paste(c(
    sprintf("NA * %s", items[[1]]),
    sprintf("%s * %s", llabels.between, items)),
    collapse = " + ")
  variances.between <- paste(sprintf(
    "%s~~%s*%s", items, rlabels.between, items),
    collapse = "\n")

  ## if only two items, need to add constraints
  ## for proper estimation
  if (identical(length(items), 2L)) {
    constraints.within <- paste0(
      constraints.within,
      "\nwl1 == wl2\n")
    constraints.between <- paste0(
      constraints.between,
      "\nbl1 == bl2\n")
  }

  model.within <- sprintf(
   "
## within level first
level: 1
 ## single factor model
 f_within =~ %s

 ## set variances
 f_within~~1*f_within
 %s
 ## set constraints
 %s

 ## define new parameters
 num_within := (%s)^2
 denom_within := (%s)^2 + (%s)
 omega_within := num_within / denom_within
",
loadings.within,
variances.within,
constraints.within,
paste(llabels.within, collapse = " + "),
paste(llabels.within, collapse = " + "),
paste(rlabels.within, collapse = " + "))

 model.between <- sprintf(
   "
## between level second
level: 2
 ## single factor model
 f_between =~ %s

 ## set variances
 f_between~~1*f_between
 %s
 ## set constraints
 %s

 ## define new parameters
 num_between := (%s)^2
 denom_between := (%s)^2 + (%s)
 omega_between := num_between / denom_between
",
loadings.between,
variances.between,
constraints.between,
paste(llabels.between, collapse = " + "),
paste(llabels.between, collapse = " + "),
paste(rlabels.between, collapse = " + "))

  model <- sprintf("%s \n%s", model.within, model.between)
  fit <- lavaan::sem(model = model, data = data,
             cluster = id)
  output <- lavaan::parameterEstimates(fit)
  label = NULL # <- palliate R CMD check
  output <- subset(output, label %in% c("omega_within", "omega_between"))
  output <- output[, c("label", "est", "ci.lower", "ci.upper")]

  if (savemodel) {
    list(
    Results = output,
    Fit = fit)
  } else {
    list(
    Results = output)
  }
}

#' modelPerformance method for merMod objects
#'
#' For pseudo R2 by cluster, the squared correlation between observed
#' and predicted values for each cluster unit is returned.  For the overall model,
#' the marginal and conditional R2 are calculated as described in the references.
#'
#' @param object A model from \code{lmer}.
#' @param ... Additional arguments, not currently used.
#' @references For estimating the marginal and conditional R-squared values,
#'   see:
#'   Nakagawa, S. and Schielzeth, H. (2013) <doi:10.1111/j.2041-210x.2012.00261.x>
#'   "A general and simple method for obtaining R2 from generalized linear mixed-effects models"
#'   and also:
#'   Johnson, P. C. (2014) <doi:10.1111/2041-210X.12225>
#'   "Extension of Nakagawa & Schielzeth's R2GLMM to random slopes models"
#' @importFrom JWileymisc modelPerformance as.modelPerformance
#' @importFrom lme4 isREML ngrps
#' @importFrom stats model.matrix model.frame cor var nobs model.frame sigma confint
#' @importFrom stats AIC BIC logLik
#' @importFrom nlme VarCorr fixef
#' @importFrom extraoperators %flipIn%
#' @method modelPerformance merMod
#' @export
#' @return a named vector with the marginal and conditional R2 values,
#'   if \code{CLUSTER = FALSE}, otherwise, a data table with the pseudo R2
#'   for each cluster unit.
#'  A list with a \code{data.table} with the following elements:
#'   \describe{
#'   \item{Model}{A character string indicating the model type, here merMod}
#'   \item{Estimator}{A character string indicating whether the model was estimated with REML or ML}
#'   \item{N_Obs}{The number of observations}
#'   \item{N_Groups}{A character string indicating the number of unique units in each grouping/clustering variable.}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{BIC}{Bayesian Information Criterion}
#'   \item{LL}{log likelihood}
#'   \item{LLDF}{log likelihood degrees of freedom}
#'   \item{Sigma}{Residual standard deviation}
#'   \item{MarginalR2}{in sample variance explained by the fixed effects}
#'   \item{ConditionalR2}{in sample variance explained by the fixed and random effects}
#'   \item{MarginalF2}{Cohen's F2 effect size R2 / (1 - R2) based off the Marginal R2}
#'   \item{ConditionalF2}{Cohen's F2 effect size R2 / (1 - R2) based off the Conditional R2}
#'   }
#' @examples
#'
#' library(JWileymisc)
#' data(aces_daily, package = "JWileymisc")
#' m1 <- lme4::lmer(PosAff ~ 1 + (1 | UserID),
#'   data = aces_daily)
#' modelPerformance(m1)
#'
#' m1 <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID),
#'   data = aces_daily)
#' modelPerformance(m1)
#'
#' rm(m1)
modelPerformance.merMod <- function(object, ...) {
  if (isGLMM(object)) stop("currently glmer() models are not supported")
  REML <- isREML(object)
  ng <- ngrps(object)
  n <- nobs(object)
  idvars <- names(ng)

  X <- model.matrix(object)
  var.fe <- var(as.vector(X %*% fixef(object))) * (n - 1) / n
  var.re <- sum(sapply(VarCorr(object)[idvars], function(Sigma) {
    xvar <- rownames(Sigma)
    xvar <- sapply(xvar, function(v) colnames(X)[colnames(X) %flipIn% v])
    Z <- X[, xvar, drop = FALSE]
    sum(diag(crossprod(Z %*% Sigma, Z))) / n
  }))
  var.e <- sigma(object)^2
  var.total <- var.fe + var.re + var.e
  R2 <- c(
    "MarginalR2" = var.fe / var.total,
    "ConditionalR2" = (var.fe + var.re) / var.total)

  LL <- logLik(object)
  LLdf <- attr(LL, "df")

  out <- data.table(
    Model = as.character("merMod"),
    Estimator = as.character(ifelse(REML, "REML", "ML")),
    N_Obs =  as.numeric(n),
    N_Groups = as.character(paste(paste0(idvars, " (", ng, ")"), collapse = "; ")),
    AIC =   as.numeric(AIC(object)),
    BIC =   as.numeric(BIC(object)),
    LL =    as.numeric(LL),
    LLDF =  as.numeric(LLdf),
    Sigma = as.numeric(sigma(object)),
    MarginalR2 =    as.numeric(R2[["MarginalR2"]]),
    ConditionalR2 = as.numeric(R2[["ConditionalR2"]]),
    MarginalF2 =    as.numeric(R2[["MarginalR2"]] / (1 - R2[["MarginalR2"]])),
    ConditionalF2 = as.numeric(R2[["ConditionalR2"]] / (1 - R2[["ConditionalR2"]]))
  )
  out <- list(out)
  attr(out, "augmentClass") <- "merMod"

  as.modelPerformance(out)
}


## clear R CMD CHECK notes
if(getRversion() >= "2.15.1")  utils::globalVariables(c("DV", "Predicted"))

#' merMod method for R2
#'
#' For pseudo R2 by cluster, the squared correlation between observed
#' and predicted values for each cluster unit is returned.  For the overall model,
#' the marginal and conditional R2 are calculated as described in the references.
#'
#' @param object A model estimated by \code{lmer}.
#' @param cluster A logical whether to calculate individual pseudo R2 values by
#'   cluster unit (if \code{TRUE}) or the marginal and conditional R2 for the
#'   overall model (if \code{FALSE}, the default).
#' @param ... Added for consistency with generic. Not currently used.
#' @return a named vector with the marginal and conditional R2 values,
#'   if \code{CLUSTER = FALSE}, otherwise, a data table with the pseudo R2
#'   for each cluster unit.
#' @references For estimating the marginal and conditional R-squared values,
#'   see:
#'   Nakagawa, S. and Schielzeth, H. (2013) <doi:10.1111/j.2041-210x.2012.00261.x>
#'   "A general and simple method for obtaining R2 from generalized linear mixed-effects models"
#'   and also:
#'   Johnson, P. C. (2014) <doi:10.1111/2041-210X.12225>
#'   "Extension of Nakagawa & Schielzeth's R2GLMM to random slopes models"
#' @keywords utils
#' @method R2 merMod
#' @export
#' @importFrom JWileymisc R2
#' @importFrom stats model.frame cor var
#' @importFrom lme4 ngrps
#' @importFrom data.table as.data.table setnames
#' @examples
#'
#' library(JWileymisc)
#' data(aces_daily, package = "JWileymisc")
#' m1 <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID),
#'   data = aces_daily)
#'
#' R2(m1)
#' R2(m1, cluster = TRUE)
#'
#' hist(R2(m1, cluster = TRUE)$R2)
#'
#' rm(m1)
R2.merMod <- function(object, cluster = FALSE, ...) {
  if (isFALSE(cluster)) {
    tmp <- modelPerformance(object)
    c(
      "MarginalR2" = tmp$Performance$MarginalR2,
      "ConditionalR2" = tmp$Performance$ConditionalR2)
  } else if (isTRUE(cluster)) {
    idvars <- names(ngrps(object))
    tmpd <- cbind(data.table(
      DV = model.frame(object)[, 1],
      Predicted = fitted(object)),
      as.data.table(
        model.frame(object)[, idvars, drop = FALSE]))

    do.call(rbind, lapply(idvars, function(n) {
      out <- tmpd[, .(
        IDVariable = n,
        R2 = cor(DV, Predicted)^2), by = get(n)]
      setnames(out, old = "get", new = "ID")
      return(out)
    }))
  }
}

## clear R CMD CHECK notes
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Estimator", "N_Groups", "N_Obs", "Model", "MarginalR2", "ConditionalR2", "Chi2", "LLDF"))


#' Compare two lmer models
#'
#' This function provides fit statistics and effect sizes for
#' model comparisons.  The models must be nested.
#'
#' @param model1 A model estimated by \code{lmer}.
#' @param model2 A model estimated by \code{lmer}.
#' @param ... Additional arguments, not currently used but included to match generic.
#' @return a data table with the fit indices for each model
#' and comparing models to each other.
#' @references For estimating the marginal and conditional R-squared values,
#'   see: Nakagawa, S. and Schielzeth, H. (2013). A general and simple method
#'   for obtaining R2 from generalized linear mixed-effects models.
#'   Methods in Ecology and Evolution, 4(2), 133-142. as well as:
#'   Johnson, P. C. (2014). Extension of Nakagawa & Schielzeth's R2GLMM to
#'   random slopes models. Methods in Ecology and Evolution, 5(9), 944-946.
#' @keywords utils
#' @method modelCompare merMod
#' @export
#' @importFrom JWileymisc modelCompare as.modelCompare
#' @importFrom stats pchisq
#' @importFrom lme4 refitML isLMM isREML fixef
#' @examples
#'
#' library(JWileymisc)
#' data(aces_daily, package = "JWileymisc")
#' m1 <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID),
#'   data = aces_daily)
#' m2 <- lme4::lmer(PosAff ~ STRESS + (1 | UserID),
#'   data = aces_daily)
#' m3 <- lme4::lmer(PosAff ~ STRESS + Female + (1 | UserID),
#'   data = aces_daily)
#'
#' modelCompare(m1, m2)
#' modelCompare(m2, m3)
#'
#' rm(m1, m2, m3)
modelCompare.merMod <- function(model1, model2, ...) {

  stopifnot(isLMM(model1))
  stopifnot(isLMM(model2))

  if (!identical(nobs(model1), nobs(model2))) {
    stop(sprintf(
      "Both models must have identical number of observations\n for a valid comparison but model1 had %d and\n model2 had %d observations.", nobs(model1), nobs(model2)))
  }

  if (!identical(isREML(model1), isREML(model2))) {
    message(
      sprintf("model1 is estimated with %s but model2 is estimated with %s.\nEstimators must match. Updating to both ML.",
              ifelse(isREML(model1), "REML", "ML"),
              ifelse(isREML(model2), "REML", "ML")))
    if (isREML(model1)) {
      model1 <- refitML(model1)
    }
    if (isREML(model2)) {
      model2 <- refitML(model2)
    }
  } else if (isREML(model1) && isREML(model2)) {
    if (!identical(sort(names(fixef(model1))), sort(names(fixef(model2))))) {
      message("When using REML, the fixed effects structure must be identical,\nbut was different. Refitting with ML.")
      model1 <- refitML(model1)
      model2 <- refitML(model2)
    }
  }

  model1sum <- summary(model1)
  model2sum <- summary(model2)
  df1 <- attr(model1sum$logLik, "df")
  df2 <- attr(model2sum$logLik, "df")

  if (identical(df1, df2)) {
    stop("One model must be nested within the other")
  } else if (df1 < df2) {
    ## do nothing
  } else if (df1 > df2) {
    df3 <- df1
    model3 <- model1
    model3sum <- model1sum

    model1 <- model2
    model1sum <- model2sum
    df1 <- df2

    model2 <- model3
    df2 <- df3
    model2sum <- model3sum

    rm(df3, model3, model3sum)
  }

  perf1 <- modelPerformance(model1)$Performance
  perf2 <- modelPerformance(model2)$Performance

  out <- rbind(perf1, perf2, perf2)
  out[3, (5:13)] <- perf2[, 5:13] - perf1[, 5:13]
  out[3, Estimator := ""]
  out[3, N_Groups := ""]
  out[3, N_Obs := 0L]
  out[, Model := c("Model 1", "Model 2", "Difference")]
  out[3, MarginalF2 := out[3, MarginalR2] / (1 - out[2, MarginalR2])]
  out[3, ConditionalF2 := out[3, ConditionalR2] / (1 - out[2, ConditionalR2])]
  out[, Chi2 := c(NA_real_, NA_real_, 2 * out[3, LL])]
  out[, P := c(NA_real_, NA_real_, pchisq(2 * out[3, LL], out[3, LLDF], lower.tail = FALSE))]

  out <- list(out)
  attr(out, "augmentClass") <- "merMod"

  as.modelCompare(out)
}


## clear R CMD CHECK notes
if(getRversion() >= "2.15.1") utils::globalVariables(c("var1", "var2", "sdcor", "Type",
                                                       "FE", "RE", "Terms", "Formula"))

#' estimate detailed results per variable and effect sizes for both fixed and random effects from lmer models
#'
#' This function extends the current \code{drop1} method for
#' \code{merMod} class objects from the lme4 package. Where
#' the default method to be able to drop both fixed and random
#' effects at once.
#'
#' At the moment, the function is aimed to \code{lmer} models
#' and has very few features for \code{glmer} or \code{nlmer}
#' models. The primary motivation was to provide a way to
#' provide an overall test of whether a variable
#' \dQuote{matters}.  In multilevel data, a variable may be
#' included in both the fixed and random effects. To provide
#' an overall test of whether it matters requires jointly testing
#' the fixed and random effects. This also is needed to provide
#' an overall effect size.
#'
#' The function works by generating a formula with one specific
#' variable or \dQuote{term} removed at all levels. A model is then
#' fit on this reduced formula and compared to the full model passed
#' in. This is a complex operation for mixed effects models for several
#' reasons. Firstly, \code{R} has no default mechanism for dropping
#' terms from both the fixed and random portions. Secondly,
#' mixed effects models do not accomodate all types of models. For example,
#' if a model includes only a random slope with no random intercept,
#' if the random slope was dropped, there would be no more random effects,
#' and at that point, \code{lmer} or \code{glmer} will not run the model.
#' It is theoretically possible to instead fit the model using
#' \code{lm} or \code{glm} but this becomes more complex for certain
#' model comparisons and calculations and is not currently implemented.
#' Marginal and conditional R2 values are calculated for each term,
#' and these are used also to calculate something akin to an
#' f-squared effect size.
#'
#' This is a new function and it is important to carefully evaluate
#' the results and check that they are accurate and that they are
#' sensible. Check accuracy by viewing the model formulae for each
#' reduced model and checking that those are indeed accurate.
#' In terms of checking whether a result is sensible or not,
#' there is a large literature on the difficulty interpretting
#' main effect tests in the presence of interactions. As it is
#' challenging to detect all interactions, especially ones that are
#' made outside of \code{R} formulae, all terms are tested. However,
#' it likely does not make sense to report results from dropping a
#' main effect but keeping the interaction term, so present
#' and interpret these with caution.
#'
#' @param object A \code{merMod} class object, the fitted result of
#'   \code{lmer}.
#' @param method A character vector indicating the types of confidence
#'   intervals to calculate. One of \dQuote{Wald}, \dQuote{profile}, or
#'   \dQuote{boot}.
#' @param control A \code{lmerControl()} results used to control how
#'   models are estimated when updating.
#' @param ... Additional arguments passed to \code{confint}
#' @importFrom JWileymisc modelTest as.na as.modelTest
#' @importFrom data.table as.data.table := setnames
#' @importFrom lme4 isGLMM isNLMM isLMM isREML nobars findbars
#' @importFrom lme4 ngrps lmer glmer refitML lmerControl glmerControl
#' @importFrom stats family formula nobs update
#' @importFrom lmerTest lsmeansLT
#' @method modelTest merMod
#' @export
#' @examples
#' ## these examples are slow to run
#' library(JWileymisc)
#' m1 <- lme4::lmer(extra ~ group + (1 | ID),
#' data = sleep, REML=FALSE)
#' modelTest(m1)
#'
#' \donttest{
#' data(aces_daily, package = "JWileymisc")
#'
#' strictControl <- lme4::lmerControl(optCtrl = list(
#'    algorithm = "NLOPT_LN_NELDERMEAD",
#'    xtol_abs = 1e-10,
#'    ftol_abs = 1e-10))
#'
#' m1 <- lme4::lmer(NegAff ~ STRESS + (1 + STRESS | UserID),
#'   data = aces_daily,
#'   control = strictControl)
#' modelTest(m1, method = "profile")
#'
#' m2 <- lme4::lmer(NegAff ~ STRESS + I(STRESS^2) + (1 + STRESS | UserID),
#'   data = aces_daily, control = strictControl)
#'
#' ## might normally use more bootstraps but keeping low for faster run
#' modelTest(m2, method = "boot", nsim = 100)
#' }
modelTest.merMod <- function(object, method = c("Wald", "profile", "boot"), control, ...) {

  if (isGLMM(object) || isNLMM(object)) {
    stop("GLMMs and NLMMs are not currently supported")
  }
  if (!isLMM(object)) {
    stop("Only LMMs fit with lmer() are currently supported")
  }
  method <- match.arg(method)

  if (missing(control)) {
    control <- lmerControl(
      optimizer = object@optinfo$optimizer,
      optCtrl = object@optinfo$control)
  }

  cis <- confint(object, method = method, oldNames = FALSE, ...)
  cis2 <- data.table(
    Term = rownames(cis),
    LL = cis[,1],
    UL = cis[,2])

  res <- as.data.table(as.data.frame(VarCorr(object)))
  res[, Term := ifelse(grp == "Residual",
                       "sigma",
                ifelse(
                  is.na(var2),
                  sprintf("sd_%s|%s", var1, grp),
                  sprintf("cor_%s.%s|%s", var2, var1, grp)))]
  res <- res[, .(Term = Term, Est = sdcor)]

  fes <- data.table(
    Term = names(fixef(object)),
    Est = as.numeric(fixef(object)))

  all <- merge(
    rbind(
      cbind(res, Type = "RE"),
      cbind(fes, Type = "FE")),
    cis2,
    by = "Term", all = TRUE)

  out.res <- all[Type == "RE"]
  out.fes <- all[Type == "FE"]

  objsum <- summary(object)

  if ("Pr(>|t|)" %in% colnames(objsum$coefficients)) {
    fe.p <- data.table(
      Term = rownames(objsum$coefficients),
      Pval = objsum$coefficients[, "Pr(>|t|)"])
  } else {
    fe.p <- data.table(
      Term = rownames(objsum$coefficients),
      Pval = (1 - pnorm(abs(objsum$coefficients[, "t value"]))) * 2)
  }

  out.fes <- merge(out.fes, fe.p, by = "Term", all = TRUE)

  ## check if linear mixed model is fit with REML
  ## and if so refit it with ML
  if (isLMM(object) && isREML(object)) {
    message(paste0(
      "Parameters and CIs are based on REML, \n",
      "but modelTests requires ML not REML fit for comparisons, \n",
      "and these are used in effect sizes. Refitting."))
  }
  object <- refitML(object)

  ngrps <- ngrps(object)
  out.misc <- modelPerformance(object)

  ## get formula
  f <- formula(object)

  ## fixed effects
  fe <- nobars(f)
  fe.terms <- terms(fe)
  fe.labs <- labels(fe.terms)
  fe.intercept <- if(identical(attr(fe.terms, "intercept"), 1L)) "1" else "0"

  ## random effects
  re <- lapply(findbars(f), deparse)

  re.group <- lapply(re, function(v) {
    gsub("(^.*)\\|(.*$)", "\\2", v)
  })

  re.terms <- lapply(re, function(v) {
    v <- gsub("(^.*)(\\|.*$)", "\\1", v)
    v <- sprintf("dv ~ %s", v)
    terms(as.formula(v))
  })
  re.labs <- lapply(re.terms, labels)
  re.intercept <- lapply(re.terms, function(x) {
    if(identical(attr(x, "intercept"), 1L)) "1" else "0"
  })

  ## all terms from fixed and random effects
  all.labs <- unique(c(fe.labs, unlist(re.labs)))
  tmp <- vector("character")
  for (i in seq_along(all.labs)) {
    tmp <- c(
      tmp,
      all.labs[match(TRUE, all.labs %flipIn% all.labs[i])])
  }
  all.labs <- unique(tmp)

  labs.levels <- data.table(
    Terms = all.labs,
    FE = as.integer(vapply(all.labs,
      function(v) any(unlist(fe.labs) %flipIn% v),
      FUN.VALUE = NA)),
    RE = as.integer(vapply(all.labs,
      function(v) any(unlist(re.labs) %flipIn% v),
      FUN.VALUE = NA)))
  labs.levels[, Type := paste0(FE, RE)]
  labs.levels <- labs.levels[,
    .(Type = if(Type == "11") c("11", "01") else Type),
    by = Terms]
  labs.levels[, FE := substr(Type, 1, 1) == "1"]
  labs.levels[, RE := substr(Type, 2, 2) == "1"]

  ## formula from reduced models, dropping one term at a time
  out.f <- unlist(lapply(seq_along(labs.levels$Terms), function(i) {
    use.fe.labs <- fe.labs[!((fe.labs %flipIn% labs.levels$Terms[i]) & labs.levels$FE[i])]
    use.re.labs <- lapply(re.labs, function(x) {
      if (length(x)) {
        x[!((x %flipIn%  labs.levels$Terms[i]) & labs.levels$RE[i])]
      } else {
        x
      }
    })

    fe.built <- sprintf("%s ~ %s%s%s",
                        as.character(f)[2],
                        fe.intercept,
                        if (length(use.fe.labs)) " + " else "",
                        paste(use.fe.labs, collapse = " + "))

    re.built <- lapply(seq_along(re), function(i) {
      if (re.intercept[[i]] == "0" && !length(use.re.labs[[i]])) {
        vector("character", 0L)
      } else {
        sprintf("(%s%s%s |%s)",
                re.intercept[[i]],
                if (length(use.re.labs[[i]])) " + " else "",
                paste(use.re.labs[[i]], collapse = " + "),
                re.group[[i]])
      }
    })

    re.built <- paste(unlist(re.built), collapse = " + ")

    if (nzchar(re.built)) {
      all.built <- paste(c(fe.built, re.built), collapse = " + ")
    } else {
      all.built <- NA_character_
    }
    return(all.built)
  }))

  labs.levels[, Formula := out.f]

  testm <- lapply(out.f, function(f) {
    if (!is.na(f)) {
      if (isLMM(object)) {
        lmer(as.formula(f), data = model.frame(object), REML = FALSE,
             control = control)
      } else if (isGLMM(object)) {
        glmer(as.formula(f), data = model.frame(object),
              family = family(object), control = control)
      }
    } else {
      NA
    }
  })

  NAtemplate <- modelPerformance(object)$Performance
  NAtemplate[, names(NAtemplate) := lapply(.SD, as.na)]
  NAtemplate[, Chi2 := NA_real_]
  NAtemplate[, P := NA_real_]

  nparfull <- attr(logLik(object), "df")
  out.tests <- do.call(rbind, lapply(seq_along(testm), function(i) {
    objreduced <- testm[[i]]
    v <- labs.levels$Terms[[i]]
    nparreduced <- attr(logLik(objreduced), "df")

    if (nparreduced >= nparfull) {
      msg <- sprintf("The full and reduced model had %d and %d parameters, respectively.\nThe reduced model should have fewer parameters.\nThis usually happens when there are interactions with categorical variables.\nFor an explanation, see:\nhttps://joshuawiley.com/multilevelTools/articles/lmer-vignette.html", nparfull, nparreduced)
      message(msg)
      tmp <- copy(NAtemplate)
    } else {
      
      if (!isTRUE(inherits(objreduced, "merMod"))) {
        tmp <- copy(NAtemplate)
      } else {
        tmp <- modelCompare(object, objreduced)$Comparison[3]
      }
    }
    tmp$Model <- v
    
    setnames(tmp, old = "Model", new = "Variable")
    return(tmp)
  }))

  out.tests <- cbind(out.tests, labs.levels[, -(1:2)])
  out.tests[, Type := factor(paste0(FE, RE),
                             levels = c("FALSETRUE", "TRUEFALSE", "TRUETRUE"),
                             labels = c("Random", "Fixed", "Fixed + Random"))]

  out <- list(
    FixedEffects = out.fes,
    RandomEffects = out.res,
    EffectSizes = out.tests,
    OverallModel = out.misc)
  attr(out, "augmentClass") <- "merMod"

  as.modelTest(out)
}

