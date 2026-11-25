# an assessment of fits prints

    Code
      print(two)
    Output
      -- Simulated and observed data -------------------------------------------------
      2 models assessed on 99 simulated datasets each (seed 1), with 20 observations
      at 10 locations.
      
      Mean
          Model Observed Median   95% interval At or above
       constant    0.926  1.011 0.301 to 1.494       58.6%
           line    0.926  1.160 0.450 to 1.644       76.8%
      
      SD
          Model Observed Median   95% interval At or above
       constant    0.883  0.965 0.639 to 1.375       60.6%
           line    0.883  0.987 0.675 to 1.492       67.7%
      
      Minimum
          Model Observed Median    95% interval At or above
       constant   -0.274 -0.819 -2.249 to 0.181       19.2%
           line   -0.274 -0.842 -1.927 to 0.260       24.2%
      
      Maximum
          Model Observed Median   95% interval At or above
       constant    2.596  2.694 1.867 to 4.020       57.6%
           line    2.596  2.884 1.970 to 4.256       73.7%
      
      i At or above is the percentage of simulated values at or above the observed
        one. Near 0% or 100% the observed value is unusual for the model, and 50% is
        not a target.
      i The observed variogram of the residuals is outside the 95% joint envelope in
        0 (constant) and 0 (line) of 6 distance classes.
      i The observed distribution is inside the 95% joint envelope.
      
      i `plot_predictive_check()` draws the statistics, the distribution and the
        variogram.

---

    Code
      print(summary(assess_fit(constant, nsim = 19, seed = 1)))
    Output
      -- Simulated and observed data -------------------------------------------------
       Statistic Observed Median    95% interval At or above
            Mean    0.926  1.006  0.438 to 1.666       57.9%
              SD    0.883  0.942  0.688 to 1.246       68.4%
         Minimum   -0.274 -0.998 -1.812 to 0.007       10.5%
         Maximum    2.596  2.660  2.036 to 3.488       52.6%
      
      i At or above is the percentage of simulated values at or above the observed
        one. Near 0% or 100% the observed value is unusual for the model, and 50% is
        not a target.
      i The observed variogram of the residuals is outside the 95% joint envelope in
        0 of 6 distance classes.
      i The observed distribution is inside the 95% joint envelope.

