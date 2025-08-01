test_that("qqplot matches expected output visually", {
  set.seed(42)
  vdiffr::expect_doppelganger(
    "qqplot example",
    qqplot(
      data.table::data.table(
        SNP = paste0("rs", 1:100),
        P = runif(100, 0, 1)
      ),
      p_cols = "P",
      p_col_labels = "p-value"
    )
  )
})

test_that("qqplot matches expected output visually with geom_point", {
  set.seed(42)
  vdiffr::expect_doppelganger(
    "qqplot example with geom_point",
    qqplot(
      data.table::data.table(
        SNP = paste0("rs", 1:100),
        P = runif(100, 0, 1)
      ),
      p_cols = "P",
      p_col_labels = "p-value",
      geom = "point"
    )
  )
})

test_that("Stratified qqplot matches expected output visually", {
  set.seed(43)
  daf <- data.frame(p1 = c(runif(100, 0, 1e-3), runif(900, 0, 1)),
                    p2 = c(runif(100, 0, 1e-3), runif(900, 0, 1))
             )

  vdiffr::expect_doppelganger(
    "stratified_qqplot example",
    stratified_qqplot(
      daf,
      p_col = "p1",
      q_col = "p2"
    )
  )
})

test_that("Stratified qqplot matches expected output visually with geom_point", {
  set.seed(43)
  daf <- data.frame(p1 = c(runif(100, 0, 1e-3), runif(900, 0, 1)),
                    p2 = c(runif(100, 0, 1e-3), runif(900, 0, 1))
             )

  vdiffr::expect_doppelganger(
    "stratified_qqplot example with geom_point",
    stratified_qqplot(
      daf,
      p_col = "p1",
      q_col = "p2",
      geom = "point"
    ),
  )
})

test_that("Stratified qqplot matches expected output visually with custom breaks", {
  set.seed(44)
  daf <- data.frame(p1 = c(runif(100, 0, 1e-3), runif(900, 0, 1)),
                    p2 = c(runif(100, 0, 1e-3), runif(900, 0, 1))
             )

  vdiffr::expect_doppelganger(
    "stratified_qqplot example with custom breaks",
    stratified_qqplot(
      daf,
      p_col = "p1",
      q_col = "p2",
      breaks = c(0, 0.5, 0.6, 1.0)
    )
  )
})
