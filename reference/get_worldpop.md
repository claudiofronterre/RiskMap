# Download gridded population counts from WorldPop

Finds and downloads WorldPop population counts for a country and year,
for the whole population or for chosen age groups and sexes, and returns
them as a single raster of people per cell.

## Usage

``` r
get_worldpop(
  country,
  year,
  ages = NULL,
  sex = c("both", "female", "male"),
  resolution = c("100m", "1km"),
  constrained = FALSE,
  un_adjusted = FALSE,
  dest_dir = tempdir(),
  messages = getOption("RiskMap.messages", TRUE)
)
```

## Arguments

- country:

  ISO3 code of the country, for example `"LBR"` for Liberia.

- year:

  Year of the population estimates.

- ages:

  `NULL` (the default) for the whole population, or a vector
  `c(from, to)` of ages in years, for example `c(20, Inf)` for adults
  aged 20 and over or `c(0, 4)` for children under five. Both ends must
  match the WorldPop age groups (0, 1 to 4, 5 to 9, and so on in
  five-year groups, with an open-ended oldest group).

- sex:

  `"both"` (the default), `"female"` or `"male"`. Only used with `ages`.

- resolution:

  `"100m"` (the default) or `"1km"`.

- constrained:

  Logical. `FALSE` (the default) uses the unconstrained estimates, which
  spread people over all land; `TRUE` uses the constrained estimates,
  which place people only in mapped settlements. See Details.

- un_adjusted:

  Logical. If `TRUE`, national totals are adjusted to match United
  Nations population estimates. Available for the whole population of
  the unconstrained 2000 to 2020 estimates.

- dest_dir:

  Directory for the downloaded files. Files already there are not
  downloaded again. Defaults to a temporary directory.

- messages:

  Logical. Report the product and the number of files. Defaults to the
  option `RiskMap.messages`, or `TRUE` if it is not set.

## Value

A `SpatRaster` with one layer, `population`, giving the number of people
per cell.

## Details

WorldPop distributes census-based population totals over a fine grid
with statistical models that use covariates such as settlement maps,
roads and night-time lights. The products differ in how people are
placed:

- Unconstrained estimates (`constrained = FALSE`) place people in every
  land cell, so sparsely populated or uninhabited areas receive small
  counts. They are available for every year from 2000 to 2020, at 100 m
  and 1 km, for the whole population and by age and sex. They suit
  population weights on grids of a few kilometres and analyses of past
  surveys.

- Constrained estimates (`constrained = TRUE`) place people only in
  cells mapped as built-up settlements. They describe where people live
  more sharply, particularly at fine resolution, but can miss small
  rural settlements. This function uses the WorldPop release covering
  2015 to 2030.

WorldPop also produces bottom-up estimates from microcensus surveys for
a few countries, which this function does not cover.

The function queries the WorldPop catalogue at
<https://hub.worldpop.org/rest/data>. The 1 km age and sex layers of the
2000 to 2020 release are not listed in the catalogue, and their
addresses are built from WorldPop's naming pattern. When `ages` covers
several age groups or both sexes, their counts are added together.

The result is in longitude and latitude. To aggregate the counts to the
cells of a prediction grid, project them with
[`terra::project()`](https://rspatial.github.io/terra/reference/project.html)
and `method = "sum"`, sum them into larger cells with
[`terra::aggregate()`](https://rspatial.github.io/terra/reference/aggregate.html),
and pass the result to
[`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md)
as `raster`.

## References

WorldPop (www.worldpop.org). Please cite the WorldPop products used, as
listed on the WorldPop website.

## See also

[`create_grid()`](https://claudiofronterre.github.io/RiskMap/reference/create_grid.md)

## Examples

``` r
# \donttest{
# Adults aged 20 and over in Liberia in 2011, at 1 km
adults <- get_worldpop("LBR", year = 2011, ages = c(20, Inf),
                       resolution = "1km")
#> ℹ WorldPop unconstrained 2000-2020 estimates at 1km, LBR 2011, aged 20 and
#>   over, in 26 files.
terra::global(adults, "sum", na.rm = TRUE)
#>                sum
#> population 1739127
# }
```
