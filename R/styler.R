## clear R CMD CHECK notes
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
           ".",
           "Term", "P", "LL", "UL",
           "ID", "value", "Est", "Pval",
           "Variable", "Type",
           "MarginalF2", "ConditionalF2"))
}

#' Format results from a linear mixed model
#'
#' @param object A list of one (or more) models estimated from lmer
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
#' data(sleepstudy, package = "lme4")
#' m1 <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
#'   data = sleepstudy)
#' m2 <- lme4::lmer(Reaction ~ Days + I(Days^2) + (1 + Days | Subject),
#'   data = sleepstudy)
#'
#' testm1 <- modelTest(m1, method = "profile")
#' testm2 <- modelTest(m2, method = "profile")
#'
#' APAStyler(testm1)
#' APAStyler(testm1,
#'   format = list(
#'     FixedEffects = "%s, %s (%s, %s)",
#'     RandomEffects = c("%s", "%s (%s, %s)"),
#'     EffectSizes = "%s, %s; %s"),
#'   pcontrol = list(digits = 3, stars = FALSE,
#'                   includeP = TRUE, includeSign = TRUE,
#'                   dropLeadingZero = TRUE))
#'
#' rm(m1, m2, testm1, testm2)
#' }
APAStyler.modelTest.merMod <- function(object,
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
    data.table(
      Term = c(
        "Model DF",
        "N (Groups)",
        "N (Observations)",
        "logLik",
        "AIC",
        "BIC",
        "Marginal R2",
        "Marginal F2",
        "Conditional R2",
        "Conditional F2"),
      Est = c(
        as.character(x$LLDF),
        as.character(x$N_Groups),
        as.character(x$N_Obs),
        formround(x$LL, digits = digits),
        formround(x$AIC, digits = digits),
        formround(x$BIC, digits = digits),
        formround(x$MarginalR2, digits = digits),
        formround(x$MarginalF2, digits = digits),
        formround(x$ConditionalR2, digits = digits),
        formround(x$ConditionalF2, digits = digits)))
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

  fe <- .formatFE(object[["FixedEffects"]], digits)
  re <- .formatRE(object[["RandomEffects"]], digits)
  misc <- .formatMISC(object[["OverallModel"]]$Performance, digits)
  ef <- .formatEFFECT(object[["EffectSizes"]], digits)

  index.fe <- data.table(
    Index = 1:nrow(fe),
    Section = 1)
  index.re <- data.table(
    Index = 1:nrow(re),
    Section = 2)
  index.misc <- data.table(
    Index = 1:nrow(misc),
    Section = 3)
  index.ef <- data.table(
    Index = 1:nrow(ef),
    Section = 4)

  out <- rbind(fe, re, misc, ef)
  attr(out, "index") <- rbind(index.fe, index.re, index.misc, index.ef)

  return(out)
}
