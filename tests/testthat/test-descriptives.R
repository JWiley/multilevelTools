test_that("meanDecompose works with one grouping factor", {
  md <- meanDecompose(mpg ~ cyl, data = mtcars)

  ## only one grouping factor so grouping means + residuals = 2
  expect_length(md, 2L)

  ## residuals are same as original data
  expect_equal(nrow(md[[length(md)]]), nrow(mtcars))

  ## cyl means are equal to unique cyl
  expect_equal(nrow(md[[1]]), length(unique(mtcars$cyl)))

  ## residuals have mean about 0
  expect_equal(mean(md[[length(md)]]$X, na.rm = TRUE),
                    0, tolerance = .05)
})

test_that("meanDecompose works with two grouping factors", {
  md <- meanDecompose(mpg ~ vs + cyl, data = mtcars)

  ## two grouping factors + residuals = 3
  expect_length(md, 3L)

  ## residuals are same as original data
  expect_identical(nrow(md[[length(md)]]), nrow(mtcars))

  ## vs means are equal to unique vs
  expect_identical(nrow(md[[1]]), length(unique(mtcars$vs)))

  ## vs by cyl means are equal to unique vs x cyl
  expect_identical(nrow(md[[2]]), length(unique(paste0(mtcars$vs, mtcars$cyl))))

  ## residuals have mean about 0
  expect_equal(mean(md[[length(md)]]$X, na.rm = TRUE),
                    0, tolerance = .05)
})


test_that("iccMixed works with one grouping factor", {
  icc <- iccMixed("mpg", "cyl", mtcars)

  expect_identical(nrow(icc), 2L)
  expect_equal(
    sum(icc$ICC, na.rm = TRUE),
    1, tolerance = .02)
})

test_that("iccMixed works with one grouping factor and data.table", {
  icc <- iccMixed("mpg", "cyl", data.table::as.data.table(mtcars))

  expect_identical(nrow(icc), 2L)
  expect_equal(
    sum(icc$ICC, na.rm = TRUE),
    1, tolerance = .02)
})

test_that("iccMixed works with two grouping factors", {
  icc <- iccMixed("mpg", c("cyl", "am"), mtcars)

  expect_identical(nrow(icc), 3L)
  expect_equal(
    sum(icc$ICC, na.rm = TRUE),
    1, tolerance = .02)
})

test_that("iccMixed works with two grouping factors and data.table", {
  icc <- iccMixed("mpg", c("cyl", "am"), data.table::as.data.table(mtcars))

  expect_identical(nrow(icc), 3L)
  expect_equal(
    sum(icc$ICC, na.rm = TRUE),
    1, tolerance = .02)
})

test_that("iccMixed works with one grouping factor and binomial family", {
  icc <- iccMixed("am", "cyl", mtcars)

  expect_identical(nrow(icc), 2L)
  expect_equal(
    sum(icc$ICC, na.rm = TRUE),
    1, tolerance = .02)
})


test_that("nEffective works with data & gaussian (default) & data.frame", {
  nk <- nEffective(dv = "mpg", id = "cyl", data = mtcars)

  expect_identical(nrow(nk), 3L)
  expect_true(all(nk$N[1:2] <= nk$N[3]))
})

test_that("nEffective works with data & binomial & data.frame", {
  nk <- nEffective(dv = "am", id = "cyl", data = mtcars, family = "binomial")

  expect_identical(nrow(nk), 3L)
  expect_true(all(nk$N[1:2] <= nk$N[3]))
})

test_that("nEffective works with data & gaussian (default) & data.table", {
  nk <- nEffective(dv = "mpg", id = "cyl", data = data.table::as.data.table(mtcars))

  expect_identical(nrow(nk), 3L)
  expect_true(all(nk$N[1:2] <= nk$N[3]))
})

test_that("nEffective works with inputs & gaussian (default)", {
  nk <- nEffective(n = 60, k = 10, icc = .6)

  expect_identical(nrow(nk), 3L)
  expect_true(all(nk$N[1:2] <= nk$N[3]))
})

test_that("nEffective works with inputs & binomial", {
  nk <- nEffective(n = 60, k = 10, icc = .6, family = "binomial")

  expect_identical(nrow(nk), 3L)
  expect_true(all(nk$N[1:2] <= nk$N[3]))
})

test_that("meanDeviations works", {
  md <- meanDeviations(1:10)

  expect_length(md, 2)

  expect_equal(md[[1]], 5.5, tolerance = .1)

  expect_equal(mean(md[[2]], na.rm = TRUE),
    0, tolerance = .1)
})


test_that("acfByID works", {
  ## example 1
  dat <- data.table::data.table(
                       x = sin(1:30),
                       time = 1:30,
                       id = 1)

  expect_error(acfByID("x", "time", "id", data = dat, lag.max = 10))

  ac <- acfByID("x", "time", "id", data = dat, lag.max = 10L)

  ## should have 10 + 1 rows for lag 0 to 10
  expect_equal(nrow(ac), 10 + 1)
  ## correlations so expect between +/- 1
  ## add .1 for tolerance
  expect_true(all(ac$AutoCorrelation >= -1.1 & ac$AutoCorrelation <= 1.1))

  ## test that it assumes all belong to same ID if not ID specified
  ac <- acfByID("x", "time", data = dat, lag.max = 10L)

  ## should have 10 + 1 rows for lag 0 to 10
  expect_equal(nrow(ac), 10 + 1)
  ## correlations so expect between +/- 1
  ## add .1 for tolerance
  expect_true(all(ac$AutoCorrelation >= -1.1 & ac$AutoCorrelation <= 1.1))


  ## test that it works with data frame
  ac <- acfByID("x", "time", "id", data = as.data.frame(dat), lag.max = 10L)

  ## should have 10 + 1 rows for lag 0 to 10
  expect_equal(nrow(ac), 10 + 1)
  ## correlations so expect between +/- 1
  ## add .1 for tolerance
  expect_true(all(ac$AutoCorrelation >= -1.1 & ac$AutoCorrelation <= 1.1))
})


test_that(".summary.ID works", {
  x <- .summary.ID(as.data.table(mtcars), var = "mpg", idvar = "cyl")
  expect_s3_class(x, "data.table")
  expect_equal(nrow(x), 3L)
  expect_equal(ncol(x), 4L)

  x <- .summary.ID(as.data.table(mtcars), var = "mpg", idvar = "cyl", robust = TRUE)
  expect_s3_class(x, "data.table")
  expect_equal(nrow(x), 3L)
  expect_equal(ncol(x), 4L)

  expect_error(.summary.ID(mtcars))
  expect_error(.summary.ID(as.data.table(mtcars), var = c("mpg", "other"), idvar = "cyl"))
  expect_error(.summary.ID(as.data.table(mtcars), var = "mpg", idvar = c("cyl", "other")))
  expect_error(.summary.ID(as.data.table(mtcars), var = "mpg", idvar = "cyl", CI = c(.95, .99)))
  expect_error(.summary.ID(as.data.table(mtcars), var = "mpg", idvar = "mpg"))
  expect_error(.summary.ID(as.data.table(mtcars), var = "test", idvar = "mpg"))
  expect_error(.summary.ID(as.data.table(mtcars), var = "mpg", idvar = "test"))
  expect_error(.summary.ID(as.data.table(mtcars), var = "mpg", idvar = "cyl", CI = c(1.5)))
  expect_error(.summary.ID(as.data.table(mtcars), var = "mpg", idvar = "cyl", CI = "test"))
})