# Compute Unique Coordinate Identifiers

This function identifies unique coordinates from a \`sf\` (simple
feature) object and assigns an identifier to each coordinate occurrence.
It returns a list containing the identifiers for each row and a vector
of unique identifiers.

## Usage

``` r
compute_ID_coords(data_sf)
```

## Arguments

- data_sf:

  An \`sf\` object containing geometrical data from which coordinates
  are extracted.

## Value

A list with the following elements:

- ID_coords:

  An integer vector where each element corresponds to a row in the
  input, indicating the index of the unique coordinate in the full set
  of unique coordinates.

- s_unique:

  An integer vector containing the unique identifiers of all distinct
  coordinates.

## Details

The function extracts the coordinate pairs from the \`sf\` object and
determines the unique coordinates. It then assigns each row in the input
data an identifier corresponding to the unique coordinate it matches.

## Author

Emanuele Giorgi <e.giorgi@lancaster.ac.uk>
