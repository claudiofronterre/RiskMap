# an assessment of fits prints

    Code
      print(two)
    Output
      -- Simulated and observed data -------------------------------------------------
      2 models assessed on 99 simulated datasets each (seed 1), with 40 observations
      at 20 locations.
      
      Mean
          Model Observed Median   95% interval At or above
       constant    0.864  0.979 0.370 to 1.770       62.6%
           line    0.864  1.077 0.468 to 1.868       78.8%
      
      SD
          Model Observed Median   95% interval At or above
       constant    0.927  0.959 0.711 to 1.278       58.6%
           line    0.927  0.976 0.725 to 1.287       62.6%
      
      Minimum
          Model Observed Median     95% interval At or above
       constant   -0.646 -1.021 -2.252 to -0.034       20.2%
           line   -0.646 -0.996  -2.231 to 0.171       24.2%
      
      Maximum
          Model Observed Median   95% interval At or above
       constant    3.341  2.957 2.112 to 4.281       29.3%
           line    3.341  3.152 2.205 to 4.394       39.4%
      
      i At or above is the percentage of simulated values that are equal to the
        observed one or above it. An observed value outside the 95% interval is
        unusual for the model, and 50% is not a target.
      i The observed variogram of the residuals is outside the 95% joint envelope in
        0 (constant) and 0 (line) of 13 distance classes.
      i The observed distribution is inside the 95% joint envelope.
      
      i `plot_predictive_check()` draws the statistics, the distribution and the
        variogram.

---

    Code
      print(summary(assess_fit(constant, nsim = 19, seed = 1)))
    Output
      -- Simulated and observed data -------------------------------------------------
       Statistic Observed Median     95% interval At or above
            Mean    0.864  1.106   0.357 to 1.546       68.4%
              SD    0.927  0.998   0.871 to 1.193       84.2%
         Minimum   -0.646 -1.297 -2.235 to -0.274       21.1%
         Maximum    3.341  3.078   2.367 to 3.875       31.6%
      
      i At or above is the percentage of simulated values that are equal to the
        observed one or above it. An observed value outside the 95% interval is
        unusual for the model, and 50% is not a target.
      i The observed variogram of the residuals is outside the 95% joint envelope in
        0 of 13 distance classes.
      i The observed distribution is inside the 95% joint envelope.

