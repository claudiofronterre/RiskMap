test_that("check_binomial functions correctly", {

  expect_no_error(check_binomial(c(0:3)))
  expect_no_error(check_binomial(c(0,1.00000001,2)))
  expect_error(check_binomial(c(-1:2)), "'y' must only consist of zero or positive integers")
  expect_error(check_binomial(c(0,1.1,2)), "'y' must only consist of zero or positive integers")
})
