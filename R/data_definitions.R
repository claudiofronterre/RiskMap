##' @title Anopheles mosquitoes in Southern Cameroon
##' @description These data contain 116 georeferenced locations on the counts of Anopheles gambiae and
##' Anopheles coluzzii in Southern Cameroon.
##' \itemize{
##'  \item web_x x-coordinate of the spatial locations.
##'  \item web_y y-coordinate of the spatial locations.
##'  \item Locality: name of the place of the sampled location.
##'  \item An.coluzzii: counts of Anopheles coluzzii.
##'  \item An.gambiae: counts of Anopheles gambiae.
##'  \item Total: total counts of Anopheles coluzzii and Anopheles gambiae.
##'  \item elevation: elevation in meters of the sampled location.
##' }
##' The coordinate reference system is 3857.
##' @docType data
##' @keywords datasets
##' @name anopheles
##' @usage data(anopheles)
##' @format A data frame with 116 rows and 7 variables
##' @source Tene Fossog, B., Ayala, D., Acevedo, P., Kengne, P., Ngomo Abeso Mebuy,
##' I., Makanga, B., et al. (2015) Habitat segregation and ecological character displacement in cryptic African malaria
##' mosquitoes. Evolutionary Applications, 8 (4), 326-345.
NULL

##' @title Female Culex pipiens abundance (collections) in the Sacramento Metropolitan Area
##'
##' @description
##' A dataset of mosquito collection events used to quantify abundance of female
##' Culex pipiens in the Sacramento Metropolitan Area (Sacramento, Placer, El Dorado
##' counties, California, USA). Each row represents a single collection event at a
##' specific location and date.
##'
##' @format A data frame with one row per collection event (for a total of 1552 collection events) and 6 variables:
##' \describe{
##'   \item{lon}{Numeric. Longitude in decimal degrees (WGS84).}
##'   \item{lat}{Numeric. Latitude in decimal degrees (WGS84).}
##'   \item{total_females}{Integer. Total number of female Cx. pipiens captured in the event.}
##'   \item{date}{Date. Exact collection date (local).}
##'   \item{trap_nights}{Integer. Number of trap-nights for the event.}
##'   \item{trap_type}{Character. Acronym for trap type used to collect mosquitoes. These include: \code{"NJLT"} (New Jersey Light Trap), \code{"GRVD"} (Gravid Trap), and \code{"MMT"} (Mosquito Magnet Trap).}
##' }
##' @details
##' Derived from the \pkg{vectorsurvR} sample datasets by:
##' \enumerate{
##' \item Filtering to female \emph{Culex pipiens} collection records within the Sacramento
##' Metropolitan Area (SMA) using point-in-polygon against the union of Sacramento, Placer,
##' and El Dorado county boundaries.
##' \item Summing counts and trap-nights per unique \emph{location}–\emph{date}–\emph{trap type}.
##' \item Coordinates are kept as WGS84 (EPSG:4326).
##' }
##'
##' This sample dataset is intended for examples, mapping, and vector-index workflows. It is
##' not official surveillance output. Dates span the same period as the \pkg{vectorsurvR}
##' sample (e.g., 2015–2021).
##'
##' @source Constructed from \code{sample_collections} in the \code{vectorsurvR} R package.
##'
##'
##' @docType data
##' @name abund_sma
##' @aliases abund_sma
##' @usage data(abund_sma)
NULL

##' @title Heavy metal biomonitoring in Galicia
##' @description This data-set relates to two studies on lead concentration
##' in moss samples, in micrograms per gram dry weight, collected in Galicia, norther Spain. The data are from two surveys, one conducted in July 2000.
##' The variables are as follows:
##' \itemize{
##'  \item x x-coordinate of the spatial locations.
##'  \item y y-coordinate of the spatial locations.
##'  \item lead number of tested people for the presence nodules.
##' }
##' The coordinate reference system of the data is 32629.
##' @docType data
##' @keywords datasets
##' @name galicia
##' @usage data(galicia)
##' @format A data frame with 195 rows and 4 variables
##' @source Diggle, P.J., Menezes, R. and Su, T.-L. (2010).
##' Geostatistical analysis under preferential sampling (with Discussion).
##' Applied Statistics, 59, 191-232.
NULL

##' West Nile virus pool tests for female *Culex pipiens* in the Sacramento Metropolitan Area
##'
##' A dataset of PCR-tested mosquito pools used to summarize infection for female
##' *Culex pipiens* in the Sacramento Metropolitan Area (Sacramento, Placer, El Dorado
##' counties, California, USA). Each row represents a tested pool at a specific
##' location and date, with an estimated pool size.
##'
##' @format A data frame with one row per tested pool (for a total of 596 pools) and 5 variables:
##' \describe{
##'   \item{lon}{Numeric. Longitude in decimal degrees (WGS84).}
##'   \item{lat}{Numeric. Latitude in decimal degrees (WGS84).}
##'   \item{est_pool_n}{Integer. Estimated number of mosquitoes in the pool (see Details).}
##'   \item{wnv_pos}{Logical. Whether the pool tested positive for West Nile virus (\code{TRUE}/\code{FALSE}).}
##'   \item{date}{Date. Pool collection date (local).}
##' }
##'
##' @details
##' Derived from the \pkg{vectorsurvR} sample datasets by:
##' \enumerate{
##' \item Filtering to female \emph{Culex pipiens} pools with WNV testing within the SMA
##' using point-in-polygon against the union of Sacramento, Placer, and El Dorado counties.
##' \item Estimating pool size \code{est_pool_n} for each pool by:
##'   \itemize{
##'     \item Summing total female counts from nearby collection points within the \emph{same week}
##'           and within a spatial radius (e.g., 2 km) to obtain \eqn{T_{\mathrm{near}}}.
##'     \item Counting nearby pools within the same week/radius to obtain \eqn{m_{\mathrm{near}}}.
##'     \item Setting \code{est_pool_n = round(max(1, T_near / m_near))}; if no nearby collections are found,
##'           a conservative fallback (default 25) is used.
##'   }
##' \item Coordinates are kept as WGS84 (EPSG:4326).
##' }
##'
##' \strong{Important:} \code{est_pool_n} is an estimate for demonstration and pooled-modelling
##' examples; it is not the recorded laboratory pool size.
##'
##' @source Constructed from \code{sample_collections} in the \code{vectorsurvR} R package.
##'
##'
##' @docType data
##' @name infect_sma
##' @aliases infect_sma
##' @usage data(infect_sma)
NULL

##' @title Simulated data-set on the Italian peninsula
##' @description This is a simulated data-set over Italy for a continuous outcome.
##' The data-set contains 10 repeated observations for each of the 200 geo-referenced locations.
##' The variables are as follows:
##' \itemize{
##'  \item x1 ordinate of the spatial locations.
##'  \item x2 abscissa of the spatial locations.
##'  \item y simulated continuous outcome.
##'  \item region the name of the region within which a given observation falls.
##'  \item province the name of the province within which a given observation falls.
##'  \item pop_dens the population density at the location of the observation.
##'  \item ID_loc an ID identifying the location to which the observation belong.
##' }
##' The coordinate reference system of the data is 32634.
##' @docType data
##' @keywords datasets
##' @name italy_sim
##' @usage data(italy_sim)
##' @format A data frame with 2000 rows and 7 variables
NULL

##' @title River-blindness in Liberia
##' @description This data-set contains counts of reported onchocerciasis (or riverblindess)
##' cases from 91 villages sampled across Liberia.
##' The variables are as follows:
##' \itemize{
##'  \item lat latitude of the of sampled villages.
##'  \item long longitude of the sampled villages.
##'  \item ntest number of tested people for the presence nodules.
##'  \item npos number of people that tested positive for the presence of nodules.
##'  \item elevation the elevation in meters of the sampled village.
##'  \item log_elevation the log-transformed elevation in meters of the sampled village.
##' }
##' @docType data
##' @keywords datasets
##' @name liberia
##' @usage data(liberia)
##' @format A data frame with 90 rows and 6 variables
##' @source Zouré, H. G. M., Noma, M., Tekle, Afework, H., Amazigo, U. V., Diggle, P. J., Giorgi, E., and Remme, J. H. F. (2014). The Geographic Distribution of Onchocerciasis in the 20 Participating Countries of the African Programme for Onchocerciasis Control: (2) Pre-Control Endemicity Levels and Estimated Number Infected.
##'  Parasites & Vectors, 7, 326.
NULL

##' @title Loa loa prevalence data from 197 village surveys
##' @description This data-set relates to a study of the prevalence of Loa loa (eyeworm) in a series of surveys undertaken in 197 villages in west Africa (Cameroon and southern Nigeria).
##' The variables are as follows:
##'
##' \itemize{
##'   \item ROW row id: 1 to 197.
##'   \item VILLCODE village id.
##'   \item LONGITUDE Longitude in degrees.
##'   \item LATITUDE Latitude in degrees.
##'   \item NO_EXAM Number of people tested.
##'   \item NO_INF Number of positive test results.
##'   \item ELEVATION Height above sea-level in metres.
##'   \item MEAN9901 Mean of all NDVI values recorded at village location, 1999-2001
##'   \item MAX9901 Maximum of all NDVI values recorded at village location, 1999-2001
##'   \item MIN9901 Minimum of all NDVI values recorded at village location, 1999-2001
##'   \item MIN9901 Minimum of all NDVI values recorded at village location, 1999-2001
##'   \item STDEV9901 standard deviation of all NDVI values recorded at village location, 1999-2001
##' }
##'
##' @references Diggle, P.J., Thomson, M.C., Christensen, O.F., Rowlingson, B., Obsomer, V., Gardon, J., Wanji, S., Takougang, I., Enyong, P., Kamgno, J., Remme, H., Boussinesq, M. and Molyneux, D.H. (2007). Spatial modelling and prediction of Loa loa risk: decision making under uncertainty. Annals of Tropical Medicine and Parasitology, 101, 499-509.
##' @docType data
##' @keywords datasets
##' @name loaloa
##' @usage data(loaloa)
##' @format A data frame with 197 rows and 11 variables
NULL

##' @title Malaria Transmission in the Western Kenyan Highlands
##' @description The dataset contains information on 82014 individuals enrolled
##' in concurrent school and community cross-sectional surveys, conducted in 46
##' school clusters in the western Kenyan highlands. Malaria was assessed by
##' rapid diagnostic test (RDT).
##'
##' The variables are as follows:
##' \itemize{
##'  \item Cluster: unique ID for each of the 46 school clusters.
##'  \item Long: longitude coordinate of the household location.
##'  \item Lat: latitude coordinate of the household location.
##'  \item RDT: binary variable indicating the outcome of the RDT:
##'  1, if positive, and 0, if negative.
##'  \item Gender: factor variable indicating the gender of the sampled individual.
##'  \item Age: age in years of the sampled individual.
##'  \item NetUse: binary variable indicating whether the sampled individual
##'  slept under a bed net the previous night: 1, if yes, 0, if no.
##'  \item MosqCntl: binary variable indicating whether the household has used some kind
##'   of mosquito control, such as sprays and coils: 1, if yes, 0, if no.
##'  \item IRS: binary variables in indicating whether there has been indoor
##'  residual spraying (IRS) in the house in the last 12 months: 1, if yes, 0,
##'  if no.
##'  \item Travel:  binary variable indicating whether the sampled individual
##'  has traveled outside the village in the last three months: 1, if yes, 0,
##'  if no.
##'  \item SES: ordinal variable indicating the socio-economic status (SES) of
##'   the household. The variables is an integer score from 1(=poor) to 5(=rich).
##'  \item District: factor variable indicating the village of the sampled
##'  individual, "Kisii Central" or "Rachuonyo".
##'  \item Survey: factor variable indicating the survey in which the
##'  participant was enrolled, "community" or "school".
##'  \item elevation: elevation, in meters, of the recorded household location
##' }
##' @docType data
##' @keywords datasets
##' @name malkenya
##' @usage data(malkenya)
##' @format A data frame with 82014 rows and 13 variables
##' @source Stevenson, J.C., Stresman, G.H., Gitonga, C.W., Gillig, J.,
##' Owaga, C., et al. (2013). Reliability of School Surveys in Estimating Geographic
##' Variation in Malaria Transmission in the Western Kenyan Highlands.
##' PLOS ONE 8(10): e77641. doi: 10.1371/journal.pone.0077641
NULL

##' @title Covariates Dataset for Malaria Prediction in Tanzania
##'
##' @description This dataset provides covariates over a 10 by 10 km regular grid covering Tanzania. It is intended to be used together with the `tz_malaria` dataset for spatial prediction of malaria prevalence.
##'
##' \describe{
##'   \item{Population}{Population density in the area (in thousands).}
##'   \item{ITN}{Percentage of households with at least one insecticide-treated net (ITN).}
##'   \item{EVI}{Enhanced Vegetation Index, indicating vegetation density.}
##'   \item{Temperature}{Average temperature in degrees Celsius.}
##'   \item{NTL}{Nighttime light intensity, indicating urbanization and infrastructure.}
##'   \item{Precipitation}{Total precipitation in millimeters.}
##'   \item{utm_x}{UTM (Universal Transverse Mercator) x-coordinate of the grid point.}
##'   \item{utm_y}{UTM (Universal Transverse Mercator) y-coordinate of the grid point.}
##' }
##'
##' @name tz_covariates
##' @docType data
##' @usage data(tz_covariates)
##' @keywords datasets
##' @format A data frame with 8740 observations of 8 variables.
##' The CRS of the UTM coordinates is 32736.
##' @source Giorgi E, Fronterrè C, Macharia PM, Alegana VA, Snow RW, Diggle PJ. 2021 Model building and assessment of the impact of covariates for disease prevalence mapping in low-resource settings: to explain and to predict. J. R. Soc. Interface 18: 20210104. \doi{10.1098/rsif.2021.0104}
NULL

##' @title Malaria Dataset from Tanzania Demographic Health Surveys 2015
##'
##' @description This dataset contains information on malaria prevalence and associated variables from the 2015 Tanzania Demographic Health Surveys. The data includes geographical, demographic, environmental, and health-related variables.
##'
##' \describe{
##'   \item{cluster.number}{ Cluster number, identifying the survey cluster.}
##'   \item{Lat}{ Latitude of the survey cluster.}
##'   \item{Long}{ Longitude of the survey cluster.}
##'   \item{MM}{ Month of the survey (in two-digit format).}
##'   \item{YY}{ Year of the survey.}
##'   \item{UpAge}{ Upper age limit of the surveyed individuals in years.}
##'   \item{LoAge}{ Lower age limit of the surveyed individuals in years.}
##'   \item{Ex}{ Number of individuals examined for malaria.}
##'   \item{Pf}{ Number of individuals tested positive for Plasmodium falciparum (malaria parasite).}
##'   \item{PfPR2.10}{ Plasmodium falciparum parasite rate in the population (aged 2-10 years).}
##'   \item{Method}{ Method used for malaria diagnosis (e.g., Rapid Diagnostic Test (RDT)).}
##'   \item{EVI}{ Enhanced Vegetation Index, indicating vegetation density.}
##'   \item{Temperature}{ Average temperature in degrees Celsius.}
##'   \item{Precipitation}{ Total precipitation in millimeters.}
##'   \item{Population}{ Population density in the area (in thousands).}
##'   \item{ITN}{ Percentage of households with at least one insecticide-treated net (ITN).}
##'   \item{NTL}{ Nighttime light intensity, indicating urbanization and infrastructure.}
##'   \item{Urban.Rural}{ Indicator of whether the area is urban ('U') or rural ('R').}
##'   \item{utm_x}{ UTM (Universal Transverse Mercator) x-coordinate of the survey cluster.}
##'   \item{utm_y}{ UTM (Universal Transverse Mercator) y-coordinate of the survey cluster.}
##' }
##' @format A data frame with 387 rows and 20 columns, containing the following variables:
##' The CRS of the UTM coordinates is 32736.
##' @name tz_malaria
##' @docType data
##' @usage data(tz_malaria)
##' @keywords datasets
##'
##' @source \href{https://dhsprogram.com}{Tanzania Demographic Health Surveys 2015},
##' Giorgi E, Fronterrè C, Macharia PM, Alegana VA, Snow RW, Diggle PJ. (2021)
##' Model building and assessment of the impact of covariates for disease prevalence
##' mapping in low-resource settings: to explain and to predict.
##' J. R. Soc. Interface 18: 20210104. \doi{10.1098/rsif.2021.0104}
NULL


##' @title Malnutrition in Ghana
##' @description This geostatistical dataset was extracted from the Demographic and Health Survey 2014 conducted in Ghana.
##' \itemize{
##'  \item lng Longitude of the sampling cluster.
##'  \item lat Latitude of the sampling cluster.
##'  \item age age in months of the child.
##'  \item sex sex of the child.
##'  \item HAZ height-for-age Z-score.
##'  \item WAZ weight-for-age Z-score
##'  \item urb binary indicator: urban area=1; rural area=0.
##'  \item etn ethnic group.
##'  \item edu level of education of the mother, which takes integer values from 1="Poorly educated" to 3="Highly educated".
##'  \item wealth wealth score of the household, which takes integer values from 1="Poor" to 3="Rich".
##' }
##' The coordinate reference system is 3857.
##' @docType data
##' @keywords datasets
##' @name malnutrition
##' @usage data(malnutrition)
##' @format A data frame with 2671 rows and 10 variables
##' @source Demographic and Health Survey, dhsprogram.com
NULL


