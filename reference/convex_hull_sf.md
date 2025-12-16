# Convex Hull of an sf Object

Computes the convex hull of an \`sf\` object, returning the boundaries
of the smallest polygon that can enclose all geometries in the input.

## Usage

``` r
convex_hull_sf(sf_object)
```

## Arguments

- sf_object:

  An \`sf\` data frame object containing geometries.

## Value

An \`sf\` object representing the convex hull of the input geometries.

## Details

The convex hull is the smallest convex polygon that encloses all points
in the input \`sf\` object. This function computes the convex hull by
first uniting all geometries in the input using \`st_union()\`, and then
applying \`st_convex_hull()\` to obtain the polygonal boundary. The
result is returned as an \`sf\` object containing the convex hull
geometry.

## See also

[`st_convex_hull`](https://r-spatial.github.io/sf/reference/geos_unary.html),
[`st_union`](https://r-spatial.github.io/sf/reference/geos_combine.html)

## Examples

``` r
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE

# Create example sf object
points <- st_sfc(st_point(c(0,0)), st_point(c(1,1)), st_point(c(2,2)), st_point(c(0,2)))
sf_points <- st_sf(geometry = points)

# Calculate the convex hull
convex_hull_result <- convex_hull_sf(sf_points)

# Plot the result
plot(sf_points, col = 'blue', pch = 19)
plot(convex_hull_result, add = TRUE, border = 'red')
```
