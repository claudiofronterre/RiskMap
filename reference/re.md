# Random Effect Model Specification

Specifies the terms for a random effect model.

## Usage

``` r
re(...)
```

## Arguments

- ...:

  Variables representing the random effects in the model.

## Value

A list of class `re.spec` containing the following elements:

- term:

  A character vector of the specified terms.

- dim:

  The number of specified terms.

- label:

  A character string representing the full call for the random effect
  model.

## Details

The function constructs a list that includes the specified terms for the
random effects. This list can be used as a specification for a random
effect model.

## Note

At least one variable must be provided as input.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
