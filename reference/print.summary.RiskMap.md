# Print Summary of RiskMap Model

Provides a `print` method for the summary of "RiskMap" objects,
detailing the model type, parameter estimates, and other relevant
statistics.

## Usage

``` r
# S3 method for class 'summary.RiskMap'
print(x, ...)
```

## Arguments

- x:

  An object of class "summary.RiskMap".

- ...:

  other parameters.

## Value

This function is used for its side effect of printing to the console. It
does not return a value.

## Details

This function prints a detailed summary of a fitted "RiskMap" model,
including:

- The type of geostatistical model (e.g., Gaussian, Binomial, Poisson).

- Confidence intervals for parameter estimates.

- Regression coefficients with their standard errors and p-values.

- Measurement error variance, if applicable.

- Spatial process parameters, including the Matern covariance
  parameters.

- Variance of the nugget effect, if applicable.

- Unstructured random effects variances, if applicable.

- Log-likelihood of the model.

- Akaike Information Criterion (AIC) for Gaussian models.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
