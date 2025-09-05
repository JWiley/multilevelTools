

## clear R CMD CHECK notes
if(getRversion() >= "2.15.1")  utils::globalVariables(c("EffectType", "OriginalOrder"))

#' Plot Diagnostics for an lmer model
#'
#' This function creates a number of diagnostic plots
#' from lmer models.
#'
#' @param x A fitted model object from \code{lmer()}.
#' @param y Included to match the generic. Not used.
#' @param plot A logical value whether or not to plot the results or
#'   simply return the graphical objects.
#' @param ask A logical whether to ask before changing plots.
#'   Only applies to interactive environments.
#' @param ncol The number of columns to use for plots.
#'   Defaults to 1.
#' @param nrow The number of rows to use for plots.
#'   Defaults to 1.
#' @param ... Included to match the generic. Not used.
#' @return a list including plots of the residuals,
#'   residuals versus fitted values, and one list for
#'   plots of all random effects
#' @importFrom JWileymisc as.na testDistribution modelDiagnostics as.modelDiagnostics
#' @importFrom extraoperators %sle%
#' @importFrom nlme VarCorr ranef
#' @importFrom data.table as.data.table :=
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 ggtitle
#' @importFrom grDevices dev.interactive devAskNewPage
#' @method plot modelDiagnostics.merMod
#' @keywords plot
#' @export
#' @examples
#'
#' library(JWileymisc)
#' sleep[1,1] <- NA
#' m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep)
#'
#' md <- modelDiagnostics(m, ev.perc = .1)
#' md$extremeValues
#'
#' data(aces_daily, package = "JWileymisc")
#' m <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID), data = aces_daily)
#'
#' md <- modelDiagnostics(m, ev.perc = .001)
#' md$extremeValues
#' plot(md$modelDiagnostics[[2]][[2]])
#' plot(md, ncol = 2, nrow = 2)
#' plot(md, ncol = 2, nrow = 3)
#'
#' rm(m, md, sleep)
plot.modelDiagnostics.merMod <- function(x, y, plot = TRUE, ask = TRUE, ncol = 1, nrow = 1, ...) {

  presid <- plot(x$residualDiagnostics, plot = FALSE)
  pranef <- lapply(seq_along(x$modelDiagnostics), function(i) {
    plot <- plot(x$modelDiagnostics[[i]][[2]], plot = FALSE)
    plot$DensityPlot + ggtitle(x$modelDiagnostics[[i]][[1]])
  })

  ## combine residual plots with random effect plots
  p <- c(presid, pranef)

  ## total plots
  n <- length(p)

  ## total plots per figure
  k <- ncol * nrow

  ## total separate figures
  figs <- ceiling(n / k)

  if (isTRUE(plot)) {
    if (isTRUE(ask) && isTRUE(dev.interactive()) && isTRUE(figs > 1)) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for (i in 1:figs) {
      print(do.call(ggarrange, c(p[((i * k - k + 1):(i * k)) %sle% n],
                                 ncol = ncol, nrow = nrow)))
    }
  }

  return(invisible(list(
    Residuals = presid,
    RandomEffects = pranef)))
}


#' Plot Diagnostics for an lme model
#'
#' This function creates a number of diagnostic plots
#' from lme models.
#'
#' @param x A fitted model object from \code{lme()}.
#' @param y Included to match the generic. Not used.
#' @param plot A logical value whether or not to plot the results or
#'   simply return the graphical objects.
#' @param ask A logical whether to ask before changing plots.
#'   Only applies to interactive environments.
#' @param ncol The number of columns to use for plots.
#'   Defaults to 1.
#' @param nrow The number of rows to use for plots.
#'   Defaults to 1.
#' @param ... Included to match the generic. Not used.
#' @return a list including plots of the residuals,
#'   residuals versus fitted values, and one list for
#'   plots of all random effects
#' @importFrom JWileymisc as.na testDistribution modelDiagnostics as.modelDiagnostics
#' @importFrom extraoperators %sle%
#' @importFrom nlme VarCorr ranef
#' @importFrom data.table as.data.table :=
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 ggtitle
#' @importFrom grDevices dev.interactive devAskNewPage
#' @method plot modelDiagnostics.lme
#' @keywords plot
#' @export
#' @examples
#'
#' library(JWileymisc)
#' sleep[1,1] <- NA
#' m <- nlme::lme(extra ~ group, data = sleep, random = ~ 1 | ID, na.action = "na.omit")
#'
#' md <- modelDiagnostics(m, ev.perc = .1)
#' md$extremeValues
#'
#' plot(md)
#'
#' data(aces_daily, package = "JWileymisc")
#' m <- nlme::lme(PosAff ~ STRESS, data = aces_daily,
#'   random = ~ 1 + STRESS | UserID, na.action = "na.omit")
#'
#' md <- modelDiagnostics(m, ev.perc = .001)
#' md$extremeValues
#' plot(md$modelDiagnostics[[2]][[2]])
#' plot(md, ncol = 2, nrow = 2)
#' plot(md, ncol = 2, nrow = 3)
#'
#' rm(m, md, sleep)
plot.modelDiagnostics.lme <- function(x, y, plot = TRUE, ask = TRUE, ncol = 1, nrow = 1, ...) {

  presid <- plot(x$residualDiagnostics, plot = FALSE)
  pranef <- lapply(seq_along(x$modelDiagnostics), function(i) {
    plot <- plot(x$modelDiagnostics[[i]][[2]], plot = FALSE)
    plot$DensityPlot + ggtitle(x$modelDiagnostics[[i]][[1]])
  })

  ## combine residual plots with random effect plots
  p <- c(presid, pranef)

  ## total plots
  n <- length(p)

  ## total plots per figure
  k <- ncol * nrow

  ## total separate figures
  figs <- ceiling(n / k)

  if (isTRUE(plot)) {
    if (isTRUE(ask) && isTRUE(dev.interactive()) && isTRUE(figs > 1)) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for (i in 1:figs) {
      print(do.call(ggarrange, c(p[((i * k - k + 1):(i * k)) %sle% n],
                                 ncol = ncol, nrow = nrow)))
    }
  }

  return(invisible(list(
    Residuals = presid,
    RandomEffects = pranef)))
}
