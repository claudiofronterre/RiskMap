# input errors name the calling function and describe the value

    Code
      fit_example(2.5, TRUE)
    Condition
      Error in `fit_example()`:
      ! `n_iter` must be a single positive integer, not 2.5.
    Code
      fit_example(c(1, 2), TRUE)
    Condition
      Error in `fit_example()`:
      ! `n_iter` must be a single positive integer, not a double vector.
    Code
      fit_example(10, "yes")
    Condition
      Error in `fit_example()`:
      ! `verbose` must be `TRUE` or `FALSE`, not "yes".
    Code
      fit_example(10, TRUE, level = 1)
    Condition
      Error in `fit_example()`:
      ! `level` must be a single finite number strictly between 0 and 1, not 1.

