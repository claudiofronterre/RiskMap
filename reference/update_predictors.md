# Update Predictors for a RiskMap Prediction Object

This function updates the predictors of a given RiskMap prediction
object. It ensures that the new predictors match the original prediction
grid and updates the relevant components of the object accordingly.

## Usage

``` r
update_predictors(object, predictors)
```

## Arguments

- object:

  A \`RiskMap.pred.re\` object, which is the output of the
  [`pred_over_grid`](https://claudiofronterre.github.io/RiskMap/reference/pred_over_grid.md)
  function.

- predictors:

  A data frame containing the new predictor values. The number of rows
  must match the prediction grid in the \`object\`.

## Value

The updated \`RiskMap.pred.re\` object.

## Details

The function performs several checks and updates:

- Ensures that \`object\` is of class \`RiskMap.pred.re\`.

- Ensures that the number of rows in \`predictors\` matches the
  prediction grid in \`object\`.

- Removes any rows with missing values in \`predictors\` and updates the
  corresponding components of the \`object\`.

- Updates the prediction locations, the predictive samples for the
  random effects, and the linear predictor.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>

Claudio Fronterre <c.fronterr@lancaster.ac.uk>
