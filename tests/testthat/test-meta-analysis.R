test_that("effect estimates are updated with one study for elementary values", {
  dat <- data.table::data.table(
    BETA = 1, SE = 0.1,
    BETA.study1 = 2, SE.study1 = 0.1
  )

  update_meta_analysis_estimates(dat, "study1")
  expect_equal(dat$BETA, 1.5)
  expect_equal(dat$SE, 0.0707, tolerance = 1e-3)
  expect_equal(dat$Z2, 450.1359, tolerance = 1e-3)
  expect_equal(dat$P, 2.03e-45)
})


test_that("effect estimates are updated for a two-study analysis", {
  dat <- data.table::data.table(
    BETA = 1, SE = 0.1,
    BETA.study1 = 1, SE.study1 = 0.1,
    BETA.study2 = 1.1, SE.study2 = 0.2
  )
  update_meta_analysis_estimates(dat, "study1")
  expect_equal(dat$BETA, 1)
  expect_equal(dat$SE, 0.0707, tolerance = 1e-3)
  expect_equal(dat$Z2, 200.0604, tolerance = 1e-3)
  expect_equal(dat$P, 2.03e-45)

  update_meta_analysis_estimates(dat, "study2")
  expect_equal(dat$BETA, 1.0111, tolerance = 1e-3)
  expect_equal(dat$SE, 0.0667, tolerance = 1e-3)
  expect_equal(dat$Z2, 230.0278, tolerance = 1e-3)
  expect_equal(dat$P, 5.88e-52)
})

test_that("effect estimates are updated for a three-study analysis", {
  dat <- data.table::data.table(
    BETA = 1, SE = 0.1,
    BETA.study1 = 1, SE.study1 = 0.1,
    BETA.study2 = 1.1, SE.study2 = 0.2
  )
  update_meta_analysis_estimates(dat, "study1")
  update_meta_analysis_estimates(dat, "study2")
  expect_equal(dat$BETA, 1.0111, tolerance = 1e-3)
  expect_equal(dat$SE, 0.0667, tolerance = 1e-3)
  expect_equal(dat$Z2, 230.0278, tolerance = 1e-3)
  expect_equal(dat$P, 5.88e-52)
})

test_that("effect estimates are updated for a three-study analysis with no initial estimates", {
  dat <- data.table::data.table(
    BETA.study1 = 1, SE.study1 = 0.1,
    BETA.study2 = 1.1, SE.study2 = 0.2,
    BETA.study3 = 1, SE.study3 = 0.1
  )
  update_meta_analysis_estimates(dat, "study1")
  update_meta_analysis_estimates(dat, "study2")
  update_meta_analysis_estimates(dat, "study3")
  expect_equal(dat$BETA, 1.0111, tolerance = 1e-3)
  expect_equal(dat$SE, 0.0667, tolerance = 1e-3)
  expect_equal(dat$Z2, 230.0278, tolerance = 1e-3)
  expect_equal(dat$P, 5.88e-52)
})

# NB: numbers from Matti Pirinnen's GWAS course: https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS9.html
test_that("matches ground truth from meta-analysis of real data from stroke GWAS compiled by Matti Pirinnen", {
  dat <- data.table::data.table(
    BETA.study1 = log(1.5), BETA.study2 = log(1.38), BETA.study3 = log(1.39),
    SE.study1 = 0.0916, SE.study2 = 0.0846, SE.study3 = 0.0967
  )

  update_meta_analysis_estimates(dat, "study1")
  update_meta_analysis_estimates(dat, "study2")
  update_meta_analysis_estimates(dat, "study3")

  expect_equal(dat$BETA, 0.3513526, tolerance = 1e-3)
  expect_equal(dat$SE, 0.05227737, tolerance = 1e-3)
  expect_equal(dat$P, 1.805663e-11, tolerance = 1e-3)
})