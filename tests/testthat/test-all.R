test_that("modelTest method works", {
  data(aces_daily, package = "JWileymisc")

  strictControl <- lme4::lmerControl(optCtrl = list(
                                       algorithm = "NLOPT_LN_NELDERMEAD",
                                       xtol_abs = 1e-10,
                                       ftol_abs = 1e-10))

  m1 <- lme4::lmer(NegAff ~ STRESS + (1 + STRESS | UserID),
                   data = aces_daily,
                   control = strictControl)
  expect_message(mt1 <- modelTest(m1))

  expect_is(mt1, c("modelTest.merMod", "modelTest"))
})

