# print.summary.RiskMap formats a gaussian model with every optional section

    Code
      print(summary_object)
    Output
      -- Linear geostatistical model -------------------------------------------------
      glgpm(formula = y ~ cov + gp(nugget = TRUE) + re(i), data = gaussian_data,
          family = "gaussian")
      The identity link.
      The limits are those of 95% confidence intervals.
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      Offset included in the linear predictor
      
      Measurement error
                             Estimate Lower limit Upper limit
      Measurement error var.      0.3         0.1         0.9
      
      Spatial Gaussian process, Matern correlation with kappa = 0.5
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget       0.10        0.02         0.5
      
      Unstructured random effects
                           Estimate Lower limit Upper limit
      i (random eff. var.)      0.4         0.2         0.8
      
      Log-likelihood  -12.3456
      AIC              38.6912

# print.summary.RiskMap reports a fixed measurement error variance

    Code
      print(summary_object)
    Output
      -- Linear geostatistical model -------------------------------------------------
      The identity link.
      The limits are those of 95% confidence intervals.
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      Offset included in the linear predictor
      
      Measurement error var. fixed at 0.1
      
      Spatial Gaussian process, Matern correlation with kappa = 0.5
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget       0.10        0.02         0.5
      
      Unstructured random effects
                           Estimate Lower limit Upper limit
      i (random eff. var.)      0.4         0.2         0.8
      
      Log-likelihood  -12.3456
      AIC              38.6912

# print.summary.RiskMap formats binomial and poisson models

    Code
      print(binomial_summary)
    Output
      -- Binomial geostatistical model -----------------------------------------------
      The canonical logit link.
      The limits are those of 95% confidence intervals.
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Spatial Gaussian process, Matern correlation with kappa = 0.5
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget effect fixed at 0
      
      Log-likelihood  -12.3456

---

    Code
      print(poisson_summary)
    Output
      -- Poisson geostatistical model ------------------------------------------------
      The canonical log link.
      The limits are those of 95% confidence intervals.
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Spatial Gaussian process, Matern correlation with kappa = 0.5
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget effect fixed at 0.25
      
      Log-likelihood  -12.3456

# print.summary.RiskMap reports the certificate and the checks to review

    Code
      print(binomial_summary)
    Output
      -- Binomial geostatistical model -----------------------------------------------
      The canonical logit link.
      The limits are those of 95% confidence intervals.
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Spatial Gaussian process, Matern correlation with kappa = 0.5
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget effect fixed at 0
      
      Log-likelihood  -12.3456
      
      Monte Carlo maximum likelihood
      v Certified. The Monte Carlo error of the estimate is within 0.059 statistical
        standard errors (an approximate 95% bound, tolerance 0.1).
      v No check needs a look. See `check_computation()` for the details.

---

    Code
      print(binomial_summary)
    Output
      -- Binomial geostatistical model -----------------------------------------------
      The canonical logit link.
      The limits are those of 95% confidence intervals.
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Spatial Gaussian process, Matern correlation with kappa = 0.5
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget effect fixed at 0
      
      Log-likelihood  -12.3456
      
      Monte Carlo maximum likelihood
      x Not certified. The Newton correction of the estimate did not converge on the
        validation draws, so the approximate 95% bound of 0.041 statistical standard
        errors may not hold. Refit with a larger `max_draws` in `set_control_mcml()`.
      ! 2 checks need a look. See `check_computation()`.

# the summary of a cross-validation gives no area or limit for binary outcomes

    Code
      print(mixed)
    Output
      -- Cross-validation scores -----------------------------------------------------
      Parameters are fixed at the estimates from all the data.
      
      Scores over 4 test observations
         Model   CRPS   SCRPS           PIT area  Limit
       model_a 0.3750 -0.6500             0.0400 0.0600
       model_b 0.3488 -0.5000 not shown (binary)       
      i Lower is better for the CRPS and the SCRPS. The CRPS is in the units of the
        response.
      i The PIT area of a correct model is below the limit with a probability of 95%.
      i For binary outcomes (model_b) the PIT area is not a check of calibration and
        is not shown, and the CRPS is the Brier score.

---

    Code
      print(binary)
    Output
      -- Cross-validation ------------------------------------------------------------
      2 models assessed on 2 test sets of 3 and 2 observations. Parameters are fixed
      at the estimates from all the data.
      
      Scores over 4 test observations
         Model   CRPS   SCRPS
       model_a 0.3488 -0.6500
       model_b 0.3488 -0.5000
      i Lower is better for the CRPS and the SCRPS. The CRPS is in the units of the
        response, which are counts.
      i For binary outcomes the PIT area is not a check of calibration and is not
        shown, and the CRPS is the Brier score.
      
      i `summary()` with `reference` compares the models with one of them, and
        `plot()` draws the calibration curves or maps the scores.

# print() of the summary of a cross-validation formats its scores

    Code
      print(scores)
    Output
      -- Cross-validation scores -----------------------------------------------------
      Parameters are fixed at the estimates from all the data.
      
      Scores over 4 test observations
         Model   CRPS   SCRPS PIT area  Limit
       model_a 0.3750 -0.6500   0.0400 0.0600
       model_b 0.3250 -0.5000   0.0500 0.0700
      i Lower is better for the CRPS and the SCRPS. The CRPS is in the units of the
        response.
      i The PIT area of a correct model is below the limit with a probability of 95%.

---

    Code
      print(compared)
    Output
      -- Cross-validation scores -----------------------------------------------------
      Parameters are fixed at the estimates from all the data.
      
      Scores over 4 test observations
         Model   CRPS  CRPS difference   SCRPS SCRPS difference PIT area  Limit
       model_a 0.3750                  -0.6500                    0.0400 0.0600
       model_b 0.3250 -0.0500 (0.1500) -0.5000 +0.1500 (0.2179)   0.0500 0.0700
      i Lower is better for the CRPS and the SCRPS. The CRPS is in the units of the
        response.
      i The PIT area of a correct model is below the limit with a probability of 95%.
      i Differences are from model_a, and negative ones are in favour of the other
        model. Their standard errors, in brackets, treat the observations as
        independent.

# print.summary.RiskMap_assess_simulation formats one spatial scale

    Code
      print(summary_object)
    Output
      -- Results of the simulation study ---------------------------------------------
      
      Mean squared error
        Model n_sim n_valid MSE_mean     MSE_sd
      1 model     2       2     0.15 0.04242641
      
      Classification by model
      Averages over the simulations by category
        Class Sensitivity Specificity PPV  NPV   CC
      1   low         0.6         0.8 0.6 0.90 0.75
      2  high         0.5         0.7 0.5 0.85 0.75
      Valid simulations by category and metric, of 2
        Class Sensitivity Specificity PPV NPV CC
      1   low           2           2   2   2  2
      2  high           2           2   2   2  2
      Proportion of correct classification across categories
        Mean 0.750, SD 0.071 and 95% simulation interval 0.702 to 0.798, from 2 of 2
        simulations.

# print.summary.RiskMap_assess_simulation labels each spatial scale

    Code
      print(summary_object)
    Output
      -- Results of the simulation study ---------------------------------------------
      
      -- Locations -------------------------------------------------------------------
      
      Mean squared error
        Model n_sim n_valid MSE_mean     MSE_sd
      1 model     2       2     0.15 0.04242641
      
      Classification by model
      Averages over the simulations by category
        Class Sensitivity Specificity PPV  NPV   CC
      1   low         0.6         0.8 0.6 0.90 0.75
      2  high         0.5         0.7 0.5 0.85 0.75
      Valid simulations by category and metric, of 2
        Class Sensitivity Specificity PPV NPV CC
      1   low           2           2   2   2  2
      2  high           2           2   2   2  2
      Proportion of correct classification across categories
        Mean 0.750, SD 0.071 and 95% simulation interval 0.702 to 0.798, from 2 of 2
        simulations.
      
      -- Areas -----------------------------------------------------------------------
      
      Mean squared error
        Model n_sim n_valid MSE_mean     MSE_sd
      1 model     2       2     0.15 0.04242641
      
      Classification by model
      Averages over the simulations by category
        Class Sensitivity Specificity PPV  NPV   CC
      1   low         0.6         0.8 0.6 0.90 0.75
      2  high         0.5         0.7 0.5 0.85 0.75
      Valid simulations by category and metric, of 2
        Class Sensitivity Specificity PPV NPV CC
      1   low           2           2   2   2  2
      2  high           2           2   2   2  2
      Proportion of correct classification across categories
        Mean 0.750, SD 0.071 and 95% simulation interval 0.702 to 0.798, from 2 of 2
        simulations.

