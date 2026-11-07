# the print of a fit without extras is unchanged

    Code
      print(fixed)
    Output
      Checks on this fit
      Fitted in 2 stages, with 2,000 fitting and 5,000 validation draws at the last
      stage.
      
      SUSPICIOUS  Variances                The nugget variance (tau2) was estimated
                                           as effectively zero, so its standard
                                           error is not reported.
                  Fix: Consider refitting without it, gp(nugget = FALSE).
      OK          Monte Carlo accuracy     The error of the estimate is within
                                           0.059 standard errors (limit 0.1).
      OK          Importance weights       The certification draws had 1,081
                                           effective draws and a tail shape of 0.46
                                           (limit 0.7).
      OK          Other maxima             The check found no higher maximum,
                                           though it cannot rule one out.
      OK          Importance distribution  The global distribution was used at
                                           every stage.
      
      These checks concern the computation, not whether the model suits the data.

