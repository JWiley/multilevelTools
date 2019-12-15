test_that("plot method for merMod model diagnostics works", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")
  m <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID), data = aces_daily)

  md <- modelDiagnostics(m, ev.perc = .01)

  pm <- plot(md, plot = FALSE)

  expect_is(pm$Residuals$ResPlot, "ggplot")
  expect_is(pm$Residuals$ResFittedPlot, "ggplot")
  expect_is(pm$RandomEffects[[1]], "ggplot")
  expect_is(pm$RandomEffects[[2]], "ggplot")
  expect_is(pm$RandomEffects[[3]], "ggplot")

  expect_invisible(plot(md, plot = TRUE, ask = FALSE, ncol = 3, nrow = 3))
})
