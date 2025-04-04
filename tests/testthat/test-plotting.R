test_that("plot method for merMod model diagnostics works", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")
  m <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID), data = aces_daily)

  md <- modelDiagnostics(m, ev.perc = .01)
  expect_s3_class(md, "modelDiagnostics.merMod")

  pm <- plot(md, plot = FALSE)

  expect_s3_class(pm$Residuals$ResPlot, "ggplot")
  expect_s3_class(pm$Residuals$ResFittedPlot, "ggplot")
  expect_s3_class(pm$RandomEffects[[1]], "ggplot")
  expect_s3_class(pm$RandomEffects[[2]], "ggplot")
  expect_s3_class(pm$RandomEffects[[3]], "ggplot")

  expect_invisible(plot(md, plot = TRUE, ask = FALSE, ncol = 3, nrow = 3))
})

test_that("plot method for merMod model diagnostics works", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")
  sleep[1,1] <- NA
  m <- nlme::lme(extra ~ group, data = sleep,
    random = ~ 1 | ID, na.action = "na.omit")

  md <- modelDiagnostics(m, ev.perc = .1)
  expect_s3_class(md, "modelDiagnostics.lme")
  expect_true(identical(nrow(md$extremeValues), 7L))

  pm <- plot(md, plot = FALSE)
  expect_s3_class(pm$Residuals$ResPlot, "ggplot")
  expect_s3_class(pm$Residuals$ResFittedPlot, "ggplot")
  expect_s3_class(pm$RandomEffects[[1]], "ggplot")

  expect_invisible(plot(md, plot = TRUE, ask = FALSE, ncol = 2, nrow = 2))
})
