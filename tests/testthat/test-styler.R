test_that("APAStyler method for merMod modelTest works", {
  skip_on_cran()

  data(aces_daily, package = "JWileymisc")
  m <- lme4::lmer(PosAff ~ STRESS + (1 + STRESS | UserID), data = aces_daily)
  expect_message(mt <- modelTest(m, method = "Wald"))

  s1 <- APAStyler(mt)

  expect_is(s1, "data.table")
  expect_identical(names(s1), c("Term", "Est"))
  expect_is(s1$Term, "character")
  expect_is(s1$Est, "character")
  expect_true(any(grepl("\\[", s1$Est)))
  expect_true(any(grepl("\\]", s1$Est)))
  expect_true(any(grepl("\\*", s1$Est)))

  s2 <- APAStyler(mt,
                  format = list(
                    FixedEffects = "%s, %s (%s, %s)",
                    RandomEffects = c("%s", "%s (%s, %s)"),
                    EffectSizes = "%s, %s; %s"),
                  pcontrol = list(digits = 3, stars = FALSE,
                                  includeP = TRUE, includeSign = TRUE,
                                  dropLeadingZero = TRUE))

  expect_is(s2, "data.table")
  expect_identical(names(s2), c("Term", "Est"))
  expect_is(s2$Term, "character")
  expect_is(s2$Est, "character")
  expect_true(any(grepl("\\(", s2$Est)))
  expect_true(any(grepl("\\)", s2$Est)))
  expect_false(any(grepl("\\*", s2$Est)))
})
