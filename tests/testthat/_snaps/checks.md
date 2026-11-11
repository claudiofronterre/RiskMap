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

# check_range leaves infinite bounds out of the message

    Code
      range_example(-5)
    Condition
      Error in `range_example()`:
      ! `max_dist` must be a single positive finite number, not -5.
    Code
      range_example(Inf)
    Condition
      Error in `range_example()`:
      ! `max_dist` must be a single positive finite number, not Inf.
    Code
      range_example(10, weight = -1)
    Condition
      Error in `range_example()`:
      ! `weight` must be a single non-negative finite number, not -1.
    Code
      range_example(10, shift = 2)
    Condition
      Error in `range_example()`:
      ! `shift` must be a single finite number less than 2, not 2.

# class and choice errors read well

    Code
      summarise_example(data.frame())
    Condition
      Error in `summarise_example()`:
      ! `object` must be an object of class <RiskMap> or <RiskMap_prediction>, not a data frame.
    Code
      summarise_example(fit, units = "miles")
    Condition
      Error in `summarise_example()`:
      ! `units` must be one of "km" or "m", not "miles".

