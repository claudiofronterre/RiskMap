# Print the Checks on a Fitted Model

Print method for objects of class `"RiskMap_check"` from
[`check_computation()`](https://claudiofronterre.github.io/RiskMap/reference/check_computation.md),
which gives the verdict. Each check has an entry, those that need a look
first. A check that needs a look is given in full, with its fix where
there is one, and a check that is `OK` takes one line, with its finding
in a few words. The extra checks that were asked for but not computed
are listed with the reason.
[`summary.RiskMap_check()`](https://claudiofronterre.github.io/RiskMap/reference/summary.RiskMap_check.md)
gives each finding in full with the numbers behind it.

## Usage

``` r
# S3 method for class 'RiskMap_check'
print(x, ...)
```

## Arguments

- x:

  An object of class `"RiskMap_check"`.

- ...:

  Not used, and an error if given, so that a misspelt argument is not
  ignored.

## Value

Invisibly returns `x`.
