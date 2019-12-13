context("residualDiagnostics for merMod")

test_that("merMod residual diagnostics fail for GLMMs", {
  skip_on_cran()

  data(cbpp, package = "lme4")

  gm1 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                     data = cbpp, family = binomial,
                     nAGQ = 0)

  expect_error(residualDiagnostics(gm1))
})

test_that("merMod residual diagnostics work for LMMs", {
  skip_on_cran()

  m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep)
  mrd <- residualDiagnostics(m)

  expect_is(mrd, c("residualDiagnostics.merMod", "residualDiagnostics"))
  expect_is(JWileymisc::as.residualDiagnostics(mrd),
            c("residualDiagnostics.merMod", "residualDiagnostics"))

  expect_equivalent(nrow(sleep), nrow(mrd$Residuals))
  expect_equivalent(
    mean(mrd$Residuals$Residuals, na.rm = TRUE),
    0,
    tolerance = .05)
})

test_that("merMod residual diagnostics work with missing data", {
  skip_on_cran()

  sleep[1,1] <- NA
  m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep,
                  na.action = na.omit)
  mrd <- residualDiagnostics(m)

  expect_is(mrd, c("residualDiagnostics.merMod", "residualDiagnostics"))
  expect_is(JWileymisc::as.residualDiagnostics(mrd),
            c("residualDiagnostics.merMod", "residualDiagnostics"))

  expect_true(nrow(mrd$Residuals) < nrow(sleep))
  expect_equivalent(
    mean(mrd$Residuals$Residuals, na.rm = TRUE),
    0,
    tolerance = .05)
})

test_that("merMod residual diagnostics work for LMMs with on the fly functions", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")
  aces_daily <- subset(aces_daily, !is.na(STRESS))
  m <- lme4::lmer(PosAff ~ poly(STRESS, 3) + (1 | UserID),
                  data = aces_daily)
  expect_warning(mrd <- residualDiagnostics(m))

  expect_is(mrd, c("residualDiagnostics.merMod", "residualDiagnostics"))
  expect_is(JWileymisc::as.residualDiagnostics(mrd),
            c("residualDiagnostics.merMod", "residualDiagnostics"))

  expect_true(nrow(mrd$Residuals) <= nrow(aces_daily))

  expect_equivalent(
    mean(mrd$Residuals$Residuals, na.rm = TRUE),
    0,
    tolerance = .05)
})
