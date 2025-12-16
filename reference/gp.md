# Gaussian Process Model Specification

Specifies the terms, smoothness, and nugget effect for a Gaussian
Process (GP) model.

## Usage

``` r
gp(..., kappa = 0.5, nugget = 0)
```

## Arguments

- ...:

  Variables representing the spatial coordinates or covariates for the
  GP model.

- kappa:

  The smoothness parameter \\\kappa\\. Default is 0.5.

- nugget:

  The nugget effect, which represents the variance of the measurement
  error. Default is 0. A positive numeric value must be provided if not
  using the default.

## Value

A list of class `gp.spec` containing the following elements:

- term:

  A character vector of the specified terms.

- kappa:

  The smoothness parameter \\\kappa\\.

- nugget:

  The nugget effect.

- dim:

  The number of specified terms.

- label:

  A character string representing the full call for the GP model.

## Details

The function constructs a list that includes the specified terms
(spatial coordinates or covariates), the smoothness parameter
\\\kappa\\, and the nugget effect. This list can be used as a
specification for a Gaussian Process model.

## Note

The nugget effect must be a positive real number if specified.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
