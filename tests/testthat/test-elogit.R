test_that("elogit matches the empirical logit formula", {
  y <- c(0, 5, 10)
  m <- c(10, 10, 10)

  expect_equal(elogit(y, m), log((y + 0.5) / (m - y + 0.5)))
  expect_true(all(is.finite(elogit(c(0, 10), c(10, 10)))))
})

test_that("elogit validates binomial inputs", {
  expect_error(elogit(11, 10), "'y' must be less than or equal to 'm'")
  expect_error(elogit(-1, 10), "'y' must contain only non-negative values")
  expect_error(elogit(1, 0), "'m' must contain only positive values")
  expect_error(elogit(c(1, 2), c(10, 10, 10)), "compatible lengths")
})
