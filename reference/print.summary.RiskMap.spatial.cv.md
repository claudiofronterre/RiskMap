# Print Summary of RiskMap Spatial Cross-Validation Scores

This function prints the matrix of cross-validation scores produced by
\`summary.RiskMap.spatial.cv\` in a readable format.

## Usage

``` r
# S3 method for class 'summary.RiskMap.spatial.cv'
print(x, ...)
```

## Arguments

- x:

  An object of class \`"summary.RiskMap.spatial.cv"\`, typically the
  output of \`summary.RiskMap.spatial.cv\`.

- ...:

  Additional arguments passed to or from other methods.

## Value

This function is used for its side effect of printing to the console. It
does not return a value.

## Details

This method is primarily used to format and display the summary score
matrix, printing it to the console. It provides a clear view of the
cross-validation performance metrics across different spatial models.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>
