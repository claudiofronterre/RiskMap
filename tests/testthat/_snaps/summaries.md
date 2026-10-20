# print.summary.RiskMap formats a gaussian model with every optional section

    Code
      print(summary_object)
    Output
      Call:
      glgpm(formula = y ~ cov + gp(nugget = TRUE) + re(i), data = gaussian_data, 
          family = "gaussian")
      
      Linear geostatistical model
      Link: identity 
      Inverse link function = x 
      
      'Lower limit' and 'Upper limit' refer to 95% confidence intervals
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      Offset included in the linear predictor
      
                             Estimate Lower limit Upper limit
      Measurement error var.      0.3         0.1         0.9
      
      Spatial Gaussian process
      Matern covariance parameters (kappa = 0.5)
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget       0.10        0.02         0.5
      
      Unstructured random effects
                           Estimate Lower limit Upper limit
      i (random eff. var.)      0.4         0.2         0.8
      
      Log-likelihood: -12.3456
      AIC: 38.6912

# print.summary.RiskMap reports a fixed measurement error variance

    Code
      print(summary_object)
    Output
      Linear geostatistical model
      Link: identity 
      Inverse link function = x 
      
      'Lower limit' and 'Upper limit' refer to 95% confidence intervals
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      Offset included in the linear predictor
      
      Measurement error var. fixed at 0.1
      
      Spatial Gaussian process
      Matern covariance parameters (kappa = 0.5)
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget       0.10        0.02         0.5
      
      Unstructured random effects
                           Estimate Lower limit Upper limit
      i (random eff. var.)      0.4         0.2         0.8
      
      Log-likelihood: -12.3456
      AIC: 38.6912

# print.summary.RiskMap formats binomial and poisson models

    Code
      print(binomial_summary)
    Output
      Binomial geostatistical model
      Link: canonical (logit) 
      Inverse link function = 1 / (1 + exp(-x)) 
      
      'Lower limit' and 'Upper limit' refer to 95% confidence intervals
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Spatial Gaussian process
      Matern covariance parameters (kappa = 0.5)
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget effect fixed at 0
      
      Log-likelihood: -12.3456

---

    Code
      print(poisson_summary)
    Output
      Poisson geostatistical model
      Link: canonical (log) 
      Inverse link function = exp(x) 
      
      'Lower limit' and 'Upper limit' refer to 95% confidence intervals
      
      Regression coefficients
                  Estimate Lower limit Upper limit StdErr z.value  p.value    
      (Intercept)     1.20        0.80        1.60   0.20     6.0 1.97e-09 ***
      cov            -0.40       -0.90        0.10   0.25    -1.6   0.1096    
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Spatial Gaussian process
      Matern covariance parameters (kappa = 0.5)
                               Estimate Lower limit Upper limit
      Spatial process var.         0.50        0.25         1.0
      Spatial corr. scale (km)     0.20        0.05         0.8
      Variance of the nugget effect fixed at 0.25
      
      Log-likelihood: -12.3456

# print.summary.RiskMap_cross_validation formats scores by model

    Code
      print(scores)
    Output
      Summary of Cross-Validation Scores
      ----------------------------------
      Model: model_a
        Test Set 1:
          CRPS: 0.2500
          SCRPS: -0.7500
        Test Set 2:
          CRPS: 0.3500
          SCRPS: -0.8500
        Overall average across test sets:
          CRPS: 0.3000
          SCRPS: -0.8000
      
      Model: model_b
        Test Set 1:
          CRPS: 0.4000
          SCRPS: -0.9000
        Test Set 2:
          CRPS: 0.5000
          SCRPS: -1.0000
        Overall average across test sets:
          CRPS: 0.4500
          SCRPS: -0.9500
      

---

    Code
      print(scores)
    Output
      Summary of Cross-Validation Scores
      ----------------------------------
      Model: model_a
        Overall average across test sets:
          CRPS: 0.3000
          SCRPS: -0.8000
      
      Model: model_b
        Overall average across test sets:
          CRPS: 0.4500
          SCRPS: -0.9500
      

# print.summary.RiskMap_assess_simulation formats one spatial scale

    Code
      print(summary_object)
    Output
      Summary of Simulation Results
      
      Mean Squared Error (MSE):
        Model n_sim n_valid MSE_mean     MSE_sd
      1 model     2       2     0.15 0.04242641
      
      Classification Results:
      
      Model: model
      
      Averages across simulations by Category:
        Class Sensitivity Specificity PPV  NPV   CC
      1   low         0.6         0.8 0.6 0.90 0.75
      2  high         0.5         0.7 0.5 0.85 0.75
      
      Number of valid simulations by Category and metric (out of 2):
        Class Sensitivity Specificity PPV NPV CC
      1   low           2           2   2   2  2
      2  high           2           2   2   2  2
      
      Proportion of Correct Classification (CC) across categories:
      Mean: 0.750, SD: 0.071, 95% simulation interval: [0.702, 0.798] (n_valid = 2 of 2)
      

# print.summary.RiskMap_assess_simulation labels each spatial scale

    Code
      print(summary_object)
    Output
      Summary of Simulation Results
      
      == Grid-level results ==
      
      Mean Squared Error (MSE):
        Model n_sim n_valid MSE_mean     MSE_sd
      1 model     2       2     0.15 0.04242641
      
      Classification Results:
      
      Model: model
      
      Averages across simulations by Category:
        Class Sensitivity Specificity PPV  NPV   CC
      1   low         0.6         0.8 0.6 0.90 0.75
      2  high         0.5         0.7 0.5 0.85 0.75
      
      Number of valid simulations by Category and metric (out of 2):
        Class Sensitivity Specificity PPV NPV CC
      1   low           2           2   2   2  2
      2  high           2           2   2   2  2
      
      Proportion of Correct Classification (CC) across categories:
      Mean: 0.750, SD: 0.071, 95% simulation interval: [0.702, 0.798] (n_valid = 2 of 2)
      
      == Area-level results ==
      
      Mean Squared Error (MSE):
        Model n_sim n_valid MSE_mean     MSE_sd
      1 model     2       2     0.15 0.04242641
      
      Classification Results:
      
      Model: model
      
      Averages across simulations by Category:
        Class Sensitivity Specificity PPV  NPV   CC
      1   low         0.6         0.8 0.6 0.90 0.75
      2  high         0.5         0.7 0.5 0.85 0.75
      
      Number of valid simulations by Category and metric (out of 2):
        Class Sensitivity Specificity PPV NPV CC
      1   low           2           2   2   2  2
      2  high           2           2   2   2  2
      
      Proportion of Correct Classification (CC) across categories:
      Mean: 0.750, SD: 0.071, 95% simulation interval: [0.702, 0.798] (n_valid = 2 of 2)
      

