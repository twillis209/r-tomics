test_that("distance_clump correctly identifies lead SNPs", {
  dat <- data.table(
    chr = c(1, 1, 2, 2, 3),
    bp = c(1000, 2000, 3000, 4000, 5000),
    p = c(1e-5, 1e-6, 1e-7, 1e-8, 1e-9)
  )

  clumped <- distance_clump(dat, "chr", "bp", "p", distance_window = 1500)

  expect_true(all(clumped$p <= 5e-8))
  expect_equal(sort(clumped$bp), c(4000, 5000))
})

test_that("distance_clump correctly identifies lead SNPs", {
  dat <- data.table(
    chr = c(1, 1, 1),
    bp = c(1000, 1100, 1501),
    p = c(1e-8, 1e-8, 1e-8)
  )

  clumped <- distance_clump(dat, "chr", "bp", "p", distance_window = 1000)

  expect_equal(sort(clumped$bp), c(1000, 1501))
})
