## clear R CMD CHECK notes
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
           "Comp", "Num", "Labels", "Ref", "Term",
           "B", "P", "LL", "UL", "SE",
           "ID", "V2", "V1",
           "value", "Est", "Pval",
           "Variable", "Type",
           "MarginalF2", "ConditionalF2"))
}

#' Format results from a linear mixed model
#'
#' @param object A list of one (or more) models estimated from lmer
#' @param modelnames An (optional) vector of names to use in
#'   the column headings for each model.
#' @param format A list giving the formatting style to be used for
#'   the fixed effecvts, random effects, and effect sizes.
#'   For the random effects, must be two options, one for when the
#'   random effects do not have confidence intervals and one when the
#'   random effects do have confidence intervals.
#' @param digits A numeric value indicating the number of digits to print.
#'   This is still in early implementation stages and currently does not
#'   change all parts of the output (which default to 2 decimals per
#'   APA style).
#' @param pcontrol A list controlling how p values are formatted.
#' @param ... Additional arguments passed to \code{confint}. Notably
#'   \code{nsim} and \code{boot.type} if the bootstrap method is used.
#' @return a data table of character data
#' @keywords misc
#' @method APAStyler modelTest.merMod
#' @export
#' @importFrom extraoperators %snin%
#' @importFrom JWileymisc formatPval APAStyler star
#' @importFrom stats pnorm
#' @importFrom nlme fixef
#' @examples
#'
#' \dontrun{
#' data(sleepstudy)
#' m1 <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
#'   data = sleepstudy)
#' m2 <- lme4::lmer(Reaction ~ Days + I(Days^2) + (1 + Days | Subject),
#'   data = sleepstudy)
#'
#' testm1 <- modelTest(m1, method = "profile")
#' testm2 <- modelTest(m2, method = "profile")
#' APAStyler(testm1)
#' APAStyler(testm1,
#'   format = list(
#'     FixedEffects = "%s, %s (%s, %s)",
#'     RandomEffects = c("%s", "%s (%s, %s)"),
#'     EffectSizes = "%s, %s; %s"),
#'   pcontrol = list(digits = 3, stars = FALSE,
#'                   includeP = TRUE, includeSign = TRUE,
#'                   dropLeadingZero = TRUE))
#' }
APAStyler.modelTest.merMod <- function(object, modelnames,
                       format = list(
                         FixedEffects = c("%s%s [%s, %s]"),
                         RandomEffects = c("%s", "%s [%s, %s]"),
                         EffectSizes = c("%s/%s, %s")),
                       digits = 2,
  pcontrol = list(digits = 3, stars = TRUE, includeP = FALSE, includeSign = FALSE,
                  dropLeadingZero = TRUE),
  ...) {

  formround <- function(x, digits) {
    format(round(x, digits = digits), digits = digits, nsmall = digits)
  }

  .formatRE <- function(x, digits) {
    tmp <- copy(x[, .(
      Term = Term,
      Est = formround(Est, digits = digits),
      LL = ifelse(is.na(LL), "", formround(LL, digits = digits)),
      UL = ifelse(is.na(UL), "", formround(UL, digits = digits)))])
    tmp[, .(
      Term = Term,
      Est = ifelse(nzchar(LL) & nzchar(UL),
                   sprintf(format$RandomEffects[2], Est, LL, UL),
                   sprintf(format$RandomEffects[1], Est)))]
    }

  .formatFE <- function(x, digits) {
    tmp <- copy(x[, .(
      Term = Term,
      Est = formround(Est, digits = digits),
      LL = ifelse(is.na(LL), "", formround(LL, digits = digits)),
      UL = ifelse(is.na(UL), "", formround(UL, digits = digits)),
      P = if (pcontrol$stars) {
            star(Pval)
          } else {
            formatPval(Pval,
                       d = pcontrol$digits,
                       sd = pcontrol$digits,
                       includeP = pcontrol$includeP,
                       includeSign = pcontrol$includeSign,
                       dropLeadingZero = pcontrol$dropLeadingZero)
          })])
    tmp[, .(
      Term = Term,
      Est = sprintf(format$FixedEffects, Est, P, LL, UL))]
  }
  .formatMISC <- function(x, digits) {
    ngrps <- gsub("^N_(.*)$", "\\1",
                  grep("^N_.*$", names(x), value = TRUE)) %snin% "Obs"
    data.table(
      Term = c(
        "Model DF",
        sprintf("N (%s)", ngrps),
        "N (Observations)",
        "logLik",
        "AIC",
        "BIC",
        "Marginal R2",
        "Conditional R2"),
      Est = c(
        as.character(x$DF),
        as.character(unlist(x[, paste0("N_", ngrps), with = FALSE])),
        as.character(x$N_Obs),
        formround(x$logLik, digits = digits),
        formround(x$AIC, digits = digits),
        formround(x$BIC, digits = digits),
        formround(x$MarginalR2, digits = digits),
        formround(x$ConditionalR2, digits = digits)))
  }

  .formatEFFECT <- function(x, digits) {
    copy(x[, .(
      Term = sprintf("%s (%s)", Variable, Type),
      Est = sprintf(format$EffectSizes,
                    formround(MarginalF2, digits = digits),
                    formround(ConditionalF2, digits = digits),
                    formatPval(P, d = pcontrol$digits, sd = pcontrol$digits,
                              includeP = TRUE, includeSign = TRUE)))])
  }



  k <- length(list)
  if (missing(modelnames)) {
    modelnames <- paste0("Model ", 1:k)
  }

  mergeIt <- function(z, func, term) {
    res <- lapply(z, function(mod) func(mod[[term]], digits = digits))
    for (i in 1:k) {
      setnames(res[[i]], old = "Est", new = modelnames[i])
    }

    if (k == 1) {
      final <- res[[1]]
    } else if (k > 1) {
      final <- merge(res[[1]], res[[2]], by = "Term", all = TRUE)
    } else if (k > 2) {
      for (i in 3:k) {
        final <- merge(final, res[[i]], by = "Term", all = TRUE)
      }
    }
    return(final)
  }

  final.fe <- mergeIt(list, .formatFE, "FixedEffects")
  final.re <- mergeIt(list, .formatRE, "RandomEffects")
  final.misc <- mergeIt(list, .formatMISC, "OverallModel")
  final.ef <- mergeIt(list, .formatEFFECT, "EffectSizes")

  placeholder <- final.fe[1]
  for (i in names(placeholder)) {
    placeholder[, (i) := ""]
  }
  place.fe <- copy(placeholder)[, Term := "Fixed Effects"]
  place.re <- copy(placeholder)[, Term := "Random Effects"]
  place.misc <- copy(placeholder)[, Term := "Overall Model"]
  place.ef <- copy(placeholder)[, Term := "Effect Sizes"]

  rbind(
    place.fe, final.fe,
    place.re, final.re,
    place.misc, final.misc,
    place.ef, final.ef)
}
