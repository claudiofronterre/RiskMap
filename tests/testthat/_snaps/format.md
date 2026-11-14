# the checks print with their symbols, in Unicode and in ASCII

    Code
      print(checks)
    Output
      -- Checks on this fit ----------------------------------------------------------
      Fitted in 2 stages, with 2,000 fitting and 5,000 validation draws at the last
      stage.
      
      ! Variances             The nugget variance (tau2) was estimated as effectively
                              zero, so its standard error is not reported.
                              > Consider refitting without it, `gp(nugget = FALSE)`.
      v Monte Carlo accuracy  Within 0.059 standard errors.
      v Importance weights    1,081 effective draws, tail shape 0.46.
      
      i `summary()` of these checks gives the evidence behind each finding.
      i These checks concern the computation, not whether the model suits the data.

---

    Code
      print(checks)
    Output
      ── Checks on this fit ──────────────────────────────────────────────────────────
      Fitted in 2 stages, with 2,000 fitting and 5,000 validation draws at the last
      stage.
      
      ! Variances             The nugget variance (tau2) was estimated as effectively
                              zero, so its standard error is not reported.
                              → Consider refitting without it, `gp(nugget = FALSE)`.
      ✔ Monte Carlo accuracy  Within 0.059 standard errors.
      ✔ Importance weights    1,081 effective draws, tail shape 0.46.
      
      ℹ `summary()` of these checks gives the evidence behind each finding.
      ℹ These checks concern the computation, not whether the model suits the data.

