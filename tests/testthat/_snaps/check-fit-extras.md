# the print gives the checks that need a look in full and the others on a line

    Code
      print(fixed)
    Output
      -- Checks on this fit ----------------------------------------------------------
      Fitted in 2 stages, with 2,000 fitting and 5,000 validation draws at the last
      stage.
      
      ! Variances                The nugget variance (tau2) was estimated as
                                 effectively zero, so its standard error is not
                                 reported.
                                 > Consider refitting without it,
                                   `gp(nugget = FALSE)`.
      v Monte Carlo accuracy     Within 0.059 standard errors.
      v Importance weights       1,081 effective draws, tail shape 0.46.
      v Other maxima             The check found no higher maximum.
      v Importance distribution  The global distribution at every stage.
      
      i `summary()` of these checks gives the evidence behind each finding.
      i These checks concern the computation, not whether the model suits the data.

