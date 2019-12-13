context("modelDiagnostics for merMod")

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

  expect_is(md, c("modelDiagnostics.merMod", "modelDiagnostics"))
  expect_is(JWileymisc::as.modelDiagnostics(md),
            c("modelDiagnostics.merMod", "modelDiagnostics"))

  expect_length(md$modelDiagnostics, 1)

  expect_equivalent(nrow(md$extremeValues), 0)
})

test_that("merMod model diagnostics work for LMMs, ev.perc = .5", {
  skip_on_cran()

  m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep)
  md <- modelDiagnostics(m, ev.perc = .5)

  expect_is(md, c("modelDiagnostics.merMod", "modelDiagnostics"))
  expect_is(JWileymisc::as.modelDiagnostics(md),
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

  expect_is(md, c("modelDiagnostics.merMod", "modelDiagnostics"))
  expect_is(JWileymisc::as.modelDiagnostics(md),
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


  m <- lme4::lmer(PosAff ~ poly(STRESS, 2) + (1 + poly(STRESS, 2) | UserID),
                  data = aces_daily, control = strictControl)

  expect_warning(md <- modelDiagnostics(m, ev.perc = .1))

  expect_is(md, c("modelDiagnostics.merMod", "modelDiagnostics"))
  expect_is(JWileymisc::as.modelDiagnostics(md),
            c("modelDiagnostics.merMod", "modelDiagnostics"))

  expect_length(md$modelDiagnostics, 4)
})
