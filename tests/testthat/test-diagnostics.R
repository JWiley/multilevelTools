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

  expect_s3_class(mrd, c("residualDiagnostics.merMod", "residualDiagnostics"))
  expect_s3_class(JWileymisc::as.residualDiagnostics(mrd),
            c("residualDiagnostics.merMod", "residualDiagnostics"))

  expect_identical(nrow(sleep), nrow(mrd$Residuals))
  expect_equal(
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

  expect_s3_class(mrd, c("residualDiagnostics.merMod", "residualDiagnostics"))
  expect_s3_class(JWileymisc::as.residualDiagnostics(mrd),
            c("residualDiagnostics.merMod", "residualDiagnostics"))

  expect_true(nrow(mrd$Residuals) < nrow(sleep))
  expect_equal(
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
  suppressWarnings(mrd <- residualDiagnostics(m))

  expect_s3_class(mrd, c("residualDiagnostics.merMod", "residualDiagnostics"))
  expect_s3_class(JWileymisc::as.residualDiagnostics(mrd),
            c("residualDiagnostics.merMod", "residualDiagnostics"))

  expect_true(nrow(mrd$Residuals) <= nrow(aces_daily))

  expect_equal(
    mean(mrd$Residuals$Residuals, na.rm = TRUE),
    0,
    tolerance = .05)
})



test_that("merMod model diagnostics fail for GLMMs", {
  skip_on_cran()

  data(cbpp, package = "lme4")

  gm1 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                     data = cbpp, family = binomial,
                     nAGQ = 0)

  expect_error(modelDiagnostics(gm1))
})

test_that("merMod model diagnostics work for LMMs, ev.perc = 0", {
  skip_on_cran()

  m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep)
  md <- modelDiagnostics(m, ev.perc = 0)

  expect_s3_class(md, c("modelDiagnostics.merMod", "modelDiagnostics"))
  expect_s3_class(JWileymisc::as.modelDiagnostics(md),
            c("modelDiagnostics.merMod", "modelDiagnostics"))

  expect_length(md$modelDiagnostics, 1)

  expect_identical(nrow(md$extremeValues), 0L)
})

test_that("merMod model diagnostics work for LMMs, ev.perc = .5", {
  skip_on_cran()

  m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep)
  md <- modelDiagnostics(m, ev.perc = .5)

  expect_s3_class(md, c("modelDiagnostics.merMod", "modelDiagnostics"))
  expect_s3_class(JWileymisc::as.modelDiagnostics(md),
            c("modelDiagnostics.merMod", "modelDiagnostics"))

  expect_length(md$modelDiagnostics, 1)

  expect_true(nrow(md$extremeValues) > 0)
})

test_that("merMod model diagnostics work with missing data", {
  skip_on_cran()

  sleep[1,1] <- NA
  m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep,
                  na.action = na.omit)
  md <- modelDiagnostics(m)

  expect_s3_class(md, c("modelDiagnostics.merMod", "modelDiagnostics"))
  expect_s3_class(JWileymisc::as.modelDiagnostics(md),
            c("modelDiagnostics.merMod", "modelDiagnostics"))

  expect_length(md$modelDiagnostics, 1)
})

test_that("merMod model diagnostics work for LMMs with on the fly functions and multiple random effects", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")
  aces_daily <- subset(aces_daily, !is.na(STRESS))

  strictControl <- lme4::lmerControl(optCtrl = list(
                                       algorithm = "NLOPT_LN_NELDERMEAD",
                                       xtol_abs = 1e-10,
                                       ftol_abs = 1e-10))


  suppressWarnings(
    m <- lme4::lmer(PosAff ~ poly(STRESS, 2) + (1 + poly(STRESS, 2) | UserID),
                  data = aces_daily, control = strictControl))

  suppressWarnings(md <- modelDiagnostics(m, ev.perc = .1))

  expect_s3_class(md, c("modelDiagnostics.merMod", "modelDiagnostics"))
  expect_s3_class(JWileymisc::as.modelDiagnostics(md),
            c("modelDiagnostics.merMod", "modelDiagnostics"))

  expect_length(md$modelDiagnostics, 4)
})
