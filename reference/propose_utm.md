# EPSG of the UTM Zone

Suggests the EPSG code for the UTM zone where the majority of the data
falls.

## Usage

``` r
propose_utm(data)
```

## Arguments

- data:

  An object of class `sf` containing the coordinates.

## Value

An integer indicating the EPSG code of the UTM zone.

## Details

The function determines the UTM zone and hemisphere where the majority
of the data points are located and proposes the corresponding EPSG code.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk> Claudio Fronterre
<c.fronterr@lancaster.ac.uk>
