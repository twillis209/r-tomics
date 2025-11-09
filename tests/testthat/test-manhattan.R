test_that("Sumstats are correctly processed for plotting", {
  dat <- data.table::data.table(
    chromosome = c(1, 1, 2, 3),
    base_pair_location = c(1, 2, 3, 4),
    p_value = c(0.1, 0.2, 0.3, 0.4)
  )

  sumstats <- process_sumstats_for_manhattan(dat)

  expect_equal(sumstats$dat$chr, c(1, 1, 2, 3))
  expect_equal(sumstats$dat$bp_cum, c(1, 2, 5 + 5e7, 10 + 1e8))
})

test_that("manhattan plot matches expected output visually", {
  set.seed(42)
  dat <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(50, 1e-10, 1e-6), runif(2150, 1e-3, 1))
  )

  sumstats <- process_sumstats_for_manhattan(dat, stat_cols = "p_value")

  vdiffr::expect_doppelganger(
    "manhattan plot example",
    draw_manhattan(
      sumstats,
      stat_col = "p_value",
      title = "Test Manhattan Plot"
    )
  )
})

test_that("manhattan plot with custom y_axis_breaks matches expected output visually", {
  set.seed(43)
  dat <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(50, 1e-10, 1e-6), runif(2150, 1e-3, 1))
  )

  sumstats <- process_sumstats_for_manhattan(dat, stat_cols = "p_value")

  vdiffr::expect_doppelganger(
    "manhattan plot with custom y_axis_breaks",
    draw_manhattan(
      sumstats,
      stat_col = "p_value",
      y_limits = c(1, 1e-10),
      y_axis_breaks = c(1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10)
    )
  )
})
