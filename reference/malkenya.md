# Malaria Transmission in the Western Kenyan Highlands

The dataset contains information on 82014 individuals enrolled in
concurrent school and community cross-sectional surveys, conducted in 46
school clusters in the western Kenyan highlands. Malaria was assessed by
rapid diagnostic test (RDT).

The variables are as follows:

- Cluster: unique ID for each of the 46 school clusters.

- Long: longitude coordinate of the household location.

- Lat: latitude coordinate of the household location.

- RDT: binary variable indicating the outcome of the RDT: 1, if
  positive, and 0, if negative.

- Gender: factor variable indicating the gender of the sampled
  individual.

- Age: age in years of the sampled individual.

- NetUse: binary variable indicating whether the sampled individual
  slept under a bed net the previous night: 1, if yes, 0, if no.

- MosqCntl: binary variable indicating whether the household has used
  some kind of mosquito control, such as sprays and coils: 1, if yes, 0,
  if no.

- IRS: binary variables in indicating whether there has been indoor
  residual spraying (IRS) in the house in the last 12 months: 1, if yes,
  0, if no.

- Travel: binary variable indicating whether the sampled individual has
  traveled outside the village in the last three months: 1, if yes, 0,
  if no.

- SES: ordinal variable indicating the socio-economic status (SES) of
  the household. The variables is an integer score from 1(=poor) to
  5(=rich).

- District: factor variable indicating the village of the sampled
  individual, "Kisii Central" or "Rachuonyo".

- Survey: factor variable indicating the survey in which the participant
  was enrolled, "community" or "school".

- elevation: elevation, in meters, of the recorded household location

## Usage

``` r
data(malkenya)
```

## Format

A data frame with 82014 rows and 13 variables

## Source

Stevenson, J.C., Stresman, G.H., Gitonga, C.W., Gillig, J., Owaga, C.,
et al. (2013). Reliability of School Surveys in Estimating Geographic
Variation in Malaria Transmission in the Western Kenyan Highlands. PLOS
ONE 8(10): e77641. doi: 10.1371/journal.pone.0077641
