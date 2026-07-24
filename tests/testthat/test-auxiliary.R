test_that("check_formula functions correctly", {

  data <- data.frame(
    y = c(1, 2, 3),
    x = c(0, 1, 2),
    z = c(0, 1, 0)
  )

  expect_no_error(check_formula(y ~ gp(x, z), data))

  expect_error(check_formula("not formula", data), "'formula' must be a 'formula'")
  expect_error(check_formula(y ~ gp(xx, z), data), "The formula term 'xx'")
  expect_error(check_formula(y ~ gp(xx, zz), data), "The formula terms 'xx', 'zz'")
  expect_error(check_formula(y ~ re(xx, zz), data), "The formula terms 'xx', 'zz'")

})
