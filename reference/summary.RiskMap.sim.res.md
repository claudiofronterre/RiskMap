# Summarize Simulation Results

Summarizes the results of model evaluations from a \`RiskMap.sim.res\`
object. Provides average metrics for classification by category and
overall correct classification (CC) summary.

## Usage

``` r
# S3 method for class 'RiskMap.sim.res'
summary(object, ...)
```

## Arguments

- object:

  An object of class \`RiskMap.sim.res\`, as returned by \`assess_sim\`.

- ...:

  Additional arguments (not used).

## Value

A list containing summary data for each model: - \`by_cat_summary\`: A
data frame with average sensitivity, specificity, PPV, NPV, and CC by
category. - \`CC_summary\`: A numeric vector with mean, 2.5th
percentile, and 97.5th percentile for CC across simulations.
