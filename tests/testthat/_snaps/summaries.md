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
      v No check needs a look. See `check_fit()` for the details.

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
        errors may not hold. Refit with a larger `max_samples` in
        `set_control_mcml()`.
      ! 2 checks need a look. See `check_fit()`.

# print.summary.RiskMap_cross_validation formats scores by model

    Code
      print(scores)
    Output
      -- Cross-validation scores -----------------------------------------------------
      
      Averages over the test sets, weighted by their sizes
         Model   CRPS   SCRPS
       model_a 0.3000 -0.8000
       model_b 0.4500 -0.9500
      
      By test set
         Model Test set   CRPS   SCRPS
       model_a        1 0.2500 -0.7500
       model_a        2 0.3500 -0.8500
       model_b        1 0.4000 -0.9000
       model_b        2 0.5000 -1.0000

---

    Code
      print(scores)
    Output
      -- Cross-validation scores -----------------------------------------------------
      
      Averages over the test sets, weighted by their sizes
         Model   CRPS   SCRPS
       model_a 0.3000 -0.8000
       model_b 0.4500 -0.9500

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
      
      -- Grid ------------------------------------------------------------------------
      
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

