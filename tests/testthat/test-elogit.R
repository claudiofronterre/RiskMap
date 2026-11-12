test_that("elogit matches the empirical logit formula", {
  y <- c(0, 5, 10)
  m <- c(10, 10, 10)

  expect_identical(elogit(y, m), log((y + 0.5) / (m - y + 0.5)))
  expect_true(all(is.finite(elogit(c(0, 10), c(10, 10)))))
})

test_that("elogit validates binomial inputs", {
  expect_input_error(elogit("a", 10), "y")
  expect_input_error(elogit(NaN, 10), "y")
  expect_input_error(elogit(11, 10), "y")
  expect_input_error(elogit(-1, 10), "y")
  expect_input_error(elogit(1, 0), "m")
  expect_error(elogit(c(1, 2), c(10, 10, 10)), "divide the longer")
})

test_that("elogit returns zero for zero length", {
  expect_identical(elogit(integer(0), double(0)), numeric(0))
})
