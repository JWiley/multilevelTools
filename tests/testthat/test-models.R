test_that("omegaSEM works", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")

  omega <- omegaSEM(
    items = c("COPEPrb", "COPEPrc", "COPEExp"),
    id = "UserID",
    data = aces_daily,
    savemodel = FALSE)

  expect_is(omega$Results, "data.frame")

  expect_true(all(
    omega$Results$est >= (omega$Results$ci.lower - .01) &
    omega$Results$est <= (omega$Results$ci.upper + .01),
    na.rm = TRUE))

  omega <- omegaSEM(
    items = c("COPEPrb", "COPEPrc", "COPEExp"),
    id = "UserID",
    data = aces_daily,
    savemodel = TRUE)

  expect_is(omega$Fit, "lavaan")
  expect_is(omega$Results, "data.frame")

  expect_true(all(
    omega$Results$est >= (omega$Results$ci.lower - .01) &
    omega$Results$est <= (omega$Results$ci.upper + .01),
    na.rm = TRUE))
})



test_that("merMod model performance fail for GLMMs", {
  skip_on_cran()

  data(cbpp, package = "lme4")

  gm1 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                     data = cbpp, family = binomial,
                     nAGQ = 0)

  expect_error(modelPerformance(gm1))
})

test_that("merMod model performance works for LMMs", {
  skip_on_cran()

  m <- lme4::lmer(extra ~ group + (1 | ID), data = sleep)
  mp <- modelPerformance(m)

  expect_is(mp, c("modelPerformance.merMod", "modelPerformance"))
  expect_is(JWileymisc::as.modelPerformance(mp),
            c("modelPerformance.merMod", "modelPerformance"))

  expect_equivalent(nrow(mp$Performance), 1)
  ## R2 should be between 0 and 1, with tolerance
  expect_true(mp$Performance$MarginalR2 >= -.05 &
              mp$Performance$MarginalR2 <= 1.05)
  expect_true(mp$Performance$ConditionalR2 >= -.05 &
              mp$Performance$ConditionalR2 <= 1.05)

  ## F2 should be equal or greater than R2, with tolerance
  expect_true(mp$Performance$MarginalF2 >=
              (mp$Performance$MarginalR2 - .01))
  expect_true(mp$Performance$ConditionalF2 >=
              (mp$Performance$ConditionalR2 - .01))
})

test_that("merMod model performance works for LMMs with intercept only", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")
  m <- lme4::lmer(PosAff ~ 1 + (1 | UserID),
                  data = aces_daily)
  mp <- modelPerformance(m)

  expect_is(mp, c("modelPerformance.merMod", "modelPerformance"))
  expect_is(JWileymisc::as.modelPerformance(mp),
            c("modelPerformance.merMod", "modelPerformance"))

  expect_equivalent(nrow(mp$Performance), 1)
  ## R2 should be between 0 and 1, with tolerance
  expect_true(mp$Performance$MarginalR2 >= -.05 &
              mp$Performance$MarginalR2 <= 1.05)
  expect_true(mp$Performance$ConditionalR2 >= -.05 &
              mp$Performance$ConditionalR2 <= 1.05)

  ## F2 should be equal or greater than R2, with tolerance
  expect_true(mp$Performance$MarginalF2 >=
              (mp$Performance$MarginalR2 - .01))
  expect_true(mp$Performance$ConditionalF2 >=
              (mp$Performance$ConditionalR2 - .01))
})

test_that("merMod model performance (& R2) works for LMMs with random slopes", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")
  m <- lme4::lmer(PosAff ~ 1 + STRESS + (1 + STRESS | UserID),
                  data = aces_daily)
  mp <- modelPerformance(m)

  expect_is(mp, c("modelPerformance.merMod", "modelPerformance"))
  expect_is(JWileymisc::as.modelPerformance(mp),
            c("modelPerformance.merMod", "modelPerformance"))

  expect_equivalent(nrow(mp$Performance), 1)
  ## R2 should be between 0 and 1, with tolerance
  expect_true(mp$Performance$MarginalR2 >= -.05 &
              mp$Performance$MarginalR2 <= 1.05)
  expect_true(mp$Performance$ConditionalR2 >= -.05 &
              mp$Performance$ConditionalR2 <= 1.05)

  ## F2 should be equal or greater than R2, with tolerance
  expect_true(mp$Performance$MarginalF2 >=
              (mp$Performance$MarginalR2 - .01))
  expect_true(mp$Performance$ConditionalF2 >=
              (mp$Performance$ConditionalR2 - .01))

  m2 <- R2(m)
  expect_is(m2, "numeric")
  ## R2 should be between 0 and 1 with tolerance
  expect_true(all(m2 >= -.05 & m2 <= 1.05, na.rm = TRUE))

  m2 <- R2(m, cluster = TRUE)
  expect_is(m2, "data.table")
})


test_that("merMod model compare works for LMMs with REML", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")

  m1 <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID),
                   data = aces_daily)
  m2 <- lme4::lmer(PosAff ~ STRESS + (1 | UserID),
                   data = aces_daily)
  m3 <- lme4::lmer(PosAff ~ STRESS + Female + (1 | UserID),
                   data = aces_daily)

  mc <- modelCompare(m1, m2)

  expect_is(mc, c("modelCompare.merMod", "modelCompare"))
  expect_is(JWileymisc::as.modelCompare(mc),
            c("modelCompare.merMod", "modelCompare"))

  expect_equivalent(nrow(mc$Comparison), 3)
  ## R2 should be between 0 and 1, with tolerance
  expect_true(all(mc$Comparison$MarginalR2 >= -.05 &
              mc$Comparison$MarginalR2 <= 1.05))
  expect_true(all(mc$Comparison$ConditionalR2 >= -.05 &
              mc$Comparison$ConditionalR2 <= 1.05))

  ## F2 should be equal or greater than R2, with tolerance
  expect_true(all(mc$Comparison$MarginalF2 >=
              (mc$Comparison$MarginalR2 - .01)))
  expect_true(all(mc$Comparison$ConditionalF2 >=
              (mc$Comparison$ConditionalR2 - .01)))


  mc <- modelCompare(m2, m3)

  expect_is(mc, c("modelCompare.merMod", "modelCompare"))
  expect_is(JWileymisc::as.modelCompare(mc),
            c("modelCompare.merMod", "modelCompare"))

  expect_equivalent(nrow(mc$Comparison), 3)
  ## R2 should be between 0 and 1, with tolerance
  expect_true(all(mc$Comparison$MarginalR2 >= -.05 &
              mc$Comparison$MarginalR2 <= 1.05))
  expect_true(all(mc$Comparison$ConditionalR2 >= -.05 &
              mc$Comparison$ConditionalR2 <= 1.05))

  ## F2 should be equal or greater than R2, with tolerance
  expect_true(all(mc$Comparison$MarginalF2 >=
              (mc$Comparison$MarginalR2 - .01)))
  expect_true(all(mc$Comparison$ConditionalF2 >=
              (mc$Comparison$ConditionalR2 - .01)))
})


test_that("merMod model test works for LMMs", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")

  strictControl <- lme4::lmerControl(optCtrl = list(
                                       algorithm = "NLOPT_LN_BOBYQA",
                                       xtol_abs = 1e-18,
                                       ftol_abs = 1e-18))

  suppressWarnings(
    m <- lme4::lmer(PosAff ~ STRESS + Female + (1 + STRESS | UserID),
                  data = aces_daily,
                  control = strictControl))

  suppressWarnings(mt <- modelTest(m))

  expect_is(mt, c("modelTest.merMod", "modelTest"))
  expect_is(JWileymisc::as.modelTest(mt),
            c("modelTest.merMod", "modelTest"))

  ## effect size R2 should be between 0 and 1, with tolerance
  expect_true(all(mt$EffectSizes$MarginalR2 >= -.05 &
              mt$EffectSizes$MarginalR2 <= 1.05))
  expect_true(all(mt$EffectSizes$ConditionalR2 >= -.05 &
              mt$EffectSizes$ConditionalR2 <= 1.05))

  ## F2 should be equal or greater than R2, with tolerance
  expect_true(all(mt$EffectSizes$MarginalF2 >=
              (mt$EffectSizes$MarginalR2 - .01)))
  expect_true(all(mt$EffectSizes$ConditionalF2 >=
              (mt$EffectSizes$ConditionalR2 - .01)))
})


