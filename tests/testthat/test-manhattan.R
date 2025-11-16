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

test_that("superimposed manhattan plot matches expected output visually", {
  set.seed(100)
  dat_before <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(50, 1e-10, 1e-6), runif(2150, 1e-3, 1))
  )

  set.seed(101)
  dat_after <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(60, 1e-12, 1e-7), runif(2140, 1e-3, 1))
  )

  sumstats_before <- process_sumstats_for_manhattan(dat_before, stat_cols = "p_value")
  sumstats_after <- process_sumstats_for_manhattan(dat_after, stat_cols = "p_value")

  vdiffr::expect_doppelganger(
    "superimposed manhattan plot",
    draw_superimposed_manhattan(
      sumstats_before,
      sumstats_after,
      stat_col = "p_value",
      alpha_before = 0.3,
      alpha_after = 1.0,
      title = "Before vs After Manhattan Plot",
      legend_labels = c("Before", "After")
    )
  )
})

test_that("superimposed manhattan plot with distinct palettes matches expected output visually", {
  set.seed(200)
  dat_before <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(50, 1e-10, 1e-6), runif(2150, 1e-3, 1))
  )

  set.seed(201)
  dat_after <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(60, 1e-12, 1e-7), runif(2140, 1e-3, 1))
  )

  sumstats_before <- process_sumstats_for_manhattan(dat_before, stat_cols = "p_value")
  sumstats_after <- process_sumstats_for_manhattan(dat_after, stat_cols = "p_value")

  vdiffr::expect_doppelganger(
    "superimposed manhattan plot with distinct palettes",
    draw_superimposed_manhattan(
      sumstats_before,
      sumstats_after,
      stat_col = "p_value",
      alpha_before = 0.4,
      alpha_after = 0.8,
      palette_before = c("#E69F00", "#F0E442"),
      palette_after = c("#0072B2", "#56B4E9"),
      title = "Before (warm) vs After (cool)",
      legend_labels = c("Before", "After")
    )
  )
})

test_that("animated manhattan plot with plotly creates valid plotly object", {
  skip_if_not_installed("plotly")
  
  set.seed(300)
  dat_before <- data.table::data.table(
    chromosome = rep(1:5, each = 50),
    base_pair_location = rep(1:50, 5) * 1e6,
    p_value = c(runif(25, 1e-8, 1e-4), runif(225, 1e-3, 1))
  )

  set.seed(301)
  dat_after <- data.table::data.table(
    chromosome = rep(1:5, each = 50),
    base_pair_location = rep(1:50, 5) * 1e6,
    p_value = c(runif(30, 1e-10, 1e-5), runif(220, 1e-3, 1))
  )

  sumstats_before <- process_sumstats_for_manhattan(
    dat_before, 
    stat_cols = "p_value",
    chromosomes = as.character(1:5)
  )
  sumstats_after <- process_sumstats_for_manhattan(
    dat_after, 
    stat_cols = "p_value",
    chromosomes = as.character(1:5)
  )

  fig <- draw_animated_manhattan_plotly(
    sumstats_before,
    sumstats_after,
    stat_col = "p_value",
    n_frames = 5,
    title = "Test Animation"
  )

  expect_s3_class(fig, "plotly")
})

test_that("animated manhattan plot with plotly and p_threshold filters correctly", {
  skip_if_not_installed("plotly")
  
  set.seed(302)
  dat_before <- data.table::data.table(
    chromosome = rep(1:5, each = 50),
    base_pair_location = rep(1:50, 5) * 1e6,
    p_value = c(runif(25, 1e-8, 1e-4), runif(225, 1e-3, 1))
  )

  set.seed(303)
  dat_after <- data.table::data.table(
    chromosome = rep(1:5, each = 50),
    base_pair_location = rep(1:50, 5) * 1e6,
    p_value = c(runif(30, 1e-10, 1e-5), runif(220, 1e-3, 1))
  )

  sumstats_before <- process_sumstats_for_manhattan(
    dat_before, 
    stat_cols = "p_value",
    chromosomes = as.character(1:5)
  )
  sumstats_after <- process_sumstats_for_manhattan(
    dat_after, 
    stat_cols = "p_value",
    chromosomes = as.character(1:5)
  )

  # Test with p_threshold
  fig <- draw_animated_manhattan_plotly(
    sumstats_before,
    sumstats_after,
    stat_col = "p_value",
    p_threshold = 0.01,
    n_frames = 5,
    title = "Test Animation with Filtering"
  )

  expect_s3_class(fig, "plotly")
})

test_that("animated manhattan plot with plotly and genome_wide_line works", {
  skip_if_not_installed("plotly")
  
  set.seed(304)
  dat_before <- data.table::data.table(
    chromosome = rep(1:5, each = 50),
    base_pair_location = rep(1:50, 5) * 1e6,
    p_value = c(runif(25, 1e-8, 1e-4), runif(225, 1e-3, 1))
  )

  set.seed(305)
  dat_after <- data.table::data.table(
    chromosome = rep(1:5, each = 50),
    base_pair_location = rep(1:50, 5) * 1e6,
    p_value = c(runif(30, 1e-10, 1e-5), runif(220, 1e-3, 1))
  )

  sumstats_before <- process_sumstats_for_manhattan(
    dat_before, 
    stat_cols = "p_value",
    chromosomes = as.character(1:5)
  )
  sumstats_after <- process_sumstats_for_manhattan(
    dat_after, 
    stat_cols = "p_value",
    chromosomes = as.character(1:5)
  )

  # Test with genome_wide_line
  fig <- draw_animated_manhattan_plotly(
    sumstats_before,
    sumstats_after,
    stat_col = "p_value",
    genome_wide_line = 5e-8,
    n_frames = 5,
    title = "Test Animation with Significance Line"
  )

  expect_s3_class(fig, "plotly")
})

test_that("back to back manhattan plot matches expected output visually", {
  set.seed(500)
  dat_top <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(50, 1e-10, 1e-6), runif(2150, 1e-3, 1))
  )

  set.seed(501)
  dat_bottom <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(60, 1e-12, 1e-7), runif(2140, 1e-3, 1))
  )

  sumstats_top <- process_sumstats_for_manhattan(dat_top, stat_cols = "p_value")
  sumstats_bottom <- process_sumstats_for_manhattan(dat_bottom, stat_cols = "p_value")

  vdiffr::expect_doppelganger(
    "back to back manhattan plot",
    draw_back_to_back_manhattan(
      sumstats_top,
      sumstats_bottom,
      stat_col = "p_value",
      title = "Back-to-Back Manhattan",
      top_label = "Dataset A",
      bottom_label = "Dataset B"
    )
  )
})

test_that("back to back manhattan plot with genome_wide_line matches expected output visually", {
  set.seed(502)
  dat_top <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(50, 1e-10, 1e-6), runif(2150, 1e-3, 1))
  )

  set.seed(503)
  dat_bottom <- data.table::data.table(
    chromosome = rep(1:22, each = 100),
    base_pair_location = rep(1:100, 22) * 1e6,
    p_value = c(runif(60, 1e-12, 1e-7), runif(2140, 1e-3, 1))
  )

  sumstats_top <- process_sumstats_for_manhattan(dat_top, stat_cols = "p_value")
  sumstats_bottom <- process_sumstats_for_manhattan(dat_bottom, stat_cols = "p_value")

  vdiffr::expect_doppelganger(
    "back to back manhattan plot with significance line",
    draw_back_to_back_manhattan(
      sumstats_top,
      sumstats_bottom,
      stat_col = "p_value",
      genome_wide_line = 5e-8,
      title = "Back-to-Back Manhattan with Threshold",
      top_label = "Dataset A",
      bottom_label = "Dataset B"
    )
  )
})

