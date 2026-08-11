##' @title Female *Culex pipiens* abundance (collections) in the Sacramento Metropolitan Area
##'
##' @description
##' Mosquito collection events used to quantify abundance of female
##' *Culex pipiens* in the Sacramento Metropolitan Area (Sacramento, Placer, El Dorado
##' counties, California, USA). Each row represents a single collection event at a
##' specific location and date.
##'
##' @format An sf object with one row per collection event (for a total of 1551 collection events) and 5 variables:
##' \describe{
##'   \item{total_females}{Integer. Total number of female Cx. pipiens captured in the event.}
##'   \item{date}{Date. Exact collection date (local).}
##'   \item{trap_nights}{Integer. Number of trap-nights for the event.}
##'   \item{trap_type}{Character. Acronym for trap type used to collect mosquitoes. These include:
##'   \code{"NJLT"} (New Jersey Light Trap), \code{"GRVD"} (Gravid Trap), and \code{"MMT"} (Mosquito Magnet Trap).}
##'   \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in degrees latitude and longitude}
##' }
##' @details
##' Derived from the \pkg{vectorsurvR} sample datasets by:
##' \enumerate{
##' \item Filtering to female *Culex pipiens* collection records within the Sacramento
##' Metropolitan Area (SMA) using point-in-polygon against the union of Sacramento, Placer,
##' and El Dorado county boundaries.
##' \item Summing counts and trap-nights per unique \emph{location}–\emph{date}–\emph{trap type}.
##' \item The coordinate reference system is 4326.
##' }
##'
##' This sample dataset is intended for examples, mapping, and vector-index workflows. It is
##' not official surveillance output. Dates span the same period as the \pkg{vectorsurvR}
##' sample (e.g., 2015–2021).
##'
##' @source Constructed from \code{sample_collections} in the \code{vectorsurvR} R package.
##'
##' @docType data
##' @name abund_sma
##' @aliases abund_sma
##' @usage data(abund_sma)
NULL


##' @title *Anopheles* mosquitoes in Southern Cameroon
##' @description Counts of *Anopheles gambiae* and
##' *Anopheles coluzzii* in Southern Cameroon.
##' @docType data
##' @keywords datasets
##' @name anopheles
##' @usage data(anopheles)
##' @format An sf object with 116 rows and 6 variables:
##' \describe{
##'  \item{Locality}{Name of the place of the sampled location.}
##'  \item{An.coluzzii}{Counts of *Anopheles coluzzii*.}
##'  \item{An.gambiae}{Counts of *Anopheles gambiae*.}
##'  \item{Total}{Total counts of *Anopheles coluzzii* and *Anopheles gambiae*.}
##'  \item{elevation}{Elevation in meters of the sampled location.}
##'  \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in degrees latitude and longitude}
##' }
##' @details The coordinate reference system is 3857.
##' @source Tene Fossog, B., Ayala, D., Acevedo, P., Kengne, P., Ngomo Abeso Mebuy,
##' I., Makanga, B., et al. (2015) Habitat segregation and ecological character displacement in cryptic African malaria
##' mosquitoes. *Evolutionary Applications*, 8 (4), 326-345.
NULL

##' @title Heavy metal biomonitoring in Galicia
##' @description Lead concentration in moss samples,
##' collected in Galicia, northern Spain. The data are from one survey conducted in July 2000.
##' @docType data
##' @keywords datasets
##' @name galicia
##' @usage data(galicia)
##' @format An sf object with 132 rows and 2 variables:
##' \describe{
##'  \item{lead}{lead concentration in micrograms per gram dry weight}
##'  \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in metres}
##' }
##' @details The coordinate reference system is 32629.
##' @source Diggle, P.J., Menezes, R. and Su, T.-L. (2010).
##' Geostatistical analysis under preferential sampling (with Discussion).
##' *Applied Statistics*, 59, 191-232.
NULL

##' @title West Nile virus pool tests for female *Culex pipiens* in the Sacramento Metropolitan Area
##'
##' @description A dataset of PCR-tested mosquito pools used to summarize infection for female
##' *Culex pipiens* in the Sacramento Metropolitan Area (Sacramento, Placer, El Dorado
##' counties, California, USA). Each row represents a tested pool at a specific
##' location and date, with an estimated pool size.
##'
##' @format An sf object with one row per tested pool (for a total of 596 pools) and 4 variables:
##' \describe{
##'   \item{est_pool_n}{Integer. Estimated number of mosquitoes in the pool (see Details).}
##'   \item{wnv_pos}{Logical. Whether the pool tested positive for West Nile virus (\code{TRUE}/\code{FALSE}).}
##'   \item{date}{Date. Pool collection date (local).}
##'   \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in degrees latitude and longitude}
##' }
##'
##' @details
##' Derived from the \pkg{vectorsurvR} sample datasets by:
##' \itemize{
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
##' \item The coordinate reference system is 4326.
##' }
##'
##' \strong{Important:} \code{est_pool_n} is an estimate for demonstration and pooled-modelling
##' examples; it is not the recorded laboratory pool size.
##'
##' @source Constructed from \code{sample_collections} in the \code{vectorsurvR} R package.
##' @docType data
##' @name infect_sma
##' @aliases infect_sma
##' @usage data(infect_sma)
NULL

##' @title Simulated dataset on the Italian peninsula
##' @description Simulated dataset over Italy for a continuous outcome.
##' @docType data
##' @keywords datasets
##' @name italy_sim
##' @usage data(italy_sim)
##' @format An sf object with 2000 rows (10 repeated observations for each of the 200 geo-referenced locations) and 6 variables:
##' \describe{
##'  \item{y}{simulated continuous outcome.}
##'  \item{region}{the name of the region within which a given observation falls.}
##'  \item{province}{the name of the province within which a given observation falls.}
##'  \item{pop_dens}{the population density at the location of the observation.}
##'  \item{ID_loc}{an ID identifying the location to which the observation belong.}
##'  \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in metres}
##' }
##' @details The coordinate reference system is 32634.
NULL

##' @title River-blindness in Liberia
##' @description Counts of reported onchocerciasis (or riverblindess)
##' cases from 90 villages sampled across Liberia.
##' @docType data
##' @keywords datasets
##' @name liberia
##' @usage data(liberia)
##' @format A data frame with 90 rows and 6 variables:
##' \describe{
##'  \item{lat}{latitude of the of sampled villages.}
##'  \item{long}{longitude of the sampled villages.}
##'  \item{ntest}{number of tested people for the presence nodules.}
##'  \item{npos}{number of people that tested positive for the presence of nodules.}
##'  \item{elevation}{the elevation in meters of the sampled village.}
##'  \item{log_elevation}{the log-transformed elevation in meters of the sampled village.}
##' }
##' @source Zouré, H. G. M., Noma, M., Tekle, Afework, H., Amazigo, U. V., Diggle,
##' P. J., Giorgi, E., and Remme, J. H. F. (2014). The Geographic Distribution of
##' Onchocerciasis in the 20 Participating Countries of the African Programme for
##' Onchocerciasis Control: (2) Pre-Control Endemicity Levels and Estimated Number Infected.
##'  *Parasites & Vectors*, 7, 326.
NULL

##' @title Loa loa prevalence data from 197 village surveys
##' @description Prevalence of Loa loa (eyeworm) in a series of surveys undertaken in
##' 197 villages in west Africa (Cameroon and southern Nigeria).
##'
##' @source Diggle, P.J., Thomson, M.C., Christensen, O.F., Rowlingson, B.,
##' Obsomer, V., Gardon, J., Wanji, S., Takougang, I., Enyong, P., Kamgno, J.,
##' Remme, H., Boussinesq, M. and Molyneux, D.H. (2007). Spatial modelling and
##' prediction of Loa loa risk: decision making under uncertainty.
##' *Annals of Tropical Medicine and Parasitology*, 101, 499-509.
##' @docType data
##' @keywords datasets
##' @name loaloa
##' @usage data(loaloa)
##' @format An sf object with 197 rows and 10 variables:
##' \describe{
##'  \item{ROW}{row id: 1 to 197.}
##'  \item{VILLCODE}{village id.}
##'  \item{NO_EXAM}{Number of people tested.}
##'  \item{NO_INF}{Number of positive test results.}
##'  \item{ELEVATION}{Height above sea-level in metres.}
##'  \item{MEAN9901}{Mean of all NDVI values recorded at village location, 1999-2001}
##'  \item{MAX9901}{Maximum of all NDVI values recorded at village location, 1999-2001}
##'  \item{MIN9901}{Minimum of all NDVI values recorded at village location, 1999-2001}
##'  \item{STDEV9901}{standard deviation of all NDVI values recorded at village location, 1999-2001}
##'  \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in degrees latitude and longitude}
##' }
##' @details The coordinate reference system is 4326.
NULL

##' @title Malaria Transmission in the Western Kenyan Highlands
##' @description The dataset contains information on 8204 individuals enrolled
##' in concurrent school and community cross-sectional surveys, conducted in 46
##' school clusters in the western Kenyan highlands. Malaria was assessed by
##' rapid diagnostic test (RDT).
##' @docType data
##' @keywords datasets
##' @name malkenya
##' @usage data(malkenya)
##' @format A data frame with 8204 rows and 13 variables:
##' \describe{
##'  \item{Cluster}{unique ID for each of the 46 school clusters.}
##'  \item{RDT}{binary variable indicating the outcome of the RDT:
##'  1, if positive, and 0, if negative.}
##'  \item{Gender}{factor variable indicating the gender of the sampled individual.}
##'  \item{Age}{age in years of the sampled individual.}
##'  \item{NetUse}{binary variable indicating whether the sampled individual
##'  slept under a bed net the previous night: 1, if yes, 0, if no.}
##'  \item{MosqCntl}{binary variable indicating whether the household has used some kind
##'   of mosquito control, such as sprays and coils: 1, if yes, 0, if no.}
##'  \item{IRS}{binary variables in indicating whether there has been indoor
##'  residual spraying (IRS) in the house in the last 12 months: 1, if yes, 0,
##'  if no.}
##'  \item{Travel}{binary variable indicating whether the sampled individual
##'  has traveled outside the village in the last three months: 1, if yes, 0,
##'  if no.}
##'  \item{SES}{ordinal variable indicating the socio-economic status (SES) of
##'   the household. The variables is an integer score from 1(=poor) to 5(=rich).}
##'  \item{District}{factor variable indicating the village of the sampled
##'  individual, "Kisii Central" or "Rachuonyo".}
##'  \item{Survey}{factor variable indicating the survey in which the
##'  participant was enrolled, "community" or "school".}
##'  \item{elevation}{elevation, in meters, of the recorded household location}
##'  \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in degrees latitude and longitude}
##' }
##' @source Stevenson, J.C., Stresman, G.H., Gitonga, C.W., Gillig, J.,
##' Owaga, C., et al. (2013). Reliability of School Surveys in Estimating Geographic
##' Variation in Malaria Transmission in the Western Kenyan Highlands.
##' *PLOS ONE* 8(10): e77641. doi: 10.1371/journal.pone.0077641
##' @details The coordinate reference system is 4326.
NULL

##' @title Malnutrition in Ghana
##' @description Malnutrition in Ghana from the Demographic and Health Survey 2014.
##' @docType data
##' @keywords datasets
##' @name malnutrition
##' @usage data(malnutrition)
##' @format A data frame with 2671 rows and 9 variables:
##' \describe{
##'  \item{age}{age in months of the child.}
##'  \item{sex}{sex of the child.}
##'  \item{HAZ}{height-for-age Z-score.}
##'  \item{WAZ}{weight-for-age Z-score.}
##'  \item{urb}{binary indicator: urban area=1; rural area=0.}
##'  \item{etn}{ethnic group.}
##'  \item{edu}{level of education of the mother, which takes integer values from 1="Poorly educated" to 3="Highly educated".}
##'  \item{wealth}{wealth score of the household, which takes integer values from 1="Poor" to 3="Rich".}
##'  \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in degrees latitude and longitude}
##' }
##' @source Demographic and Health Survey, dhsprogram.com
##' @details The coordinate reference system is 3857.
NULL


##' @title Covariates Dataset for Malaria Prediction in Tanzania
##'
##' @description This dataset provides covariates over a 10 by 10 km regular grid covering Tanzania. It is intended to be used together with the `tz_malaria` dataset for spatial prediction of malaria prevalence.
##'

##'
##' @name tz_covariates
##' @docType data
##' @usage data(tz_covariates)
##' @keywords datasets
##' @format A data frame with 8740 observations of 7 variables:
##' \describe{
##'   \item{Population}{Population density in the area (in thousands).}
##'   \item{ITN}{Percentage of households with at least one insecticide-treated net (ITN).}
##'   \item{EVI}{Enhanced Vegetation Index, indicating vegetation density.}
##'   \item{Temperature}{Average temperature in degrees Celsius.}
##'   \item{NTL}{Nighttime light intensity, indicating urbanization and infrastructure.}
##'   \item{Precipitation}{Total precipitation in millimeters.}
##'   \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in metres}
##' }
##' @details The coordinate reference system is 32736.
##' @source Giorgi E, Fronterrè C, Macharia PM, Alegana VA, Snow RW, Diggle PJ. (2021)
##'  Model building and assessment of the impact of covariates for disease prevalence
##'  mapping in low-resource settings: to explain and to predict. *J. R. Soc. Interface* 18: 20210104. \doi{10.1098/rsif.2021.0104}
NULL

##' @title Malaria Dataset from Tanzania Demographic Health Surveys 2015
##'
##' @description Malaria prevalence and associated variables from the 2015 Tanzania Demographic Health Surveys.
##' The data includes geographical, demographic, environmental, and health-related variables.
##' @format A data frame with 387 rows and 18 columns, containing the following variables:
##' \describe{
##'   \item{cluster.number}{Cluster number, identifying the survey cluster.}
##'   \item{MM}{Month of the survey (in two-digit format).}
##'   \item{YY}{Year of the survey.}
##'   \item{UpAge}{Upper age limit of the surveyed individuals in years.}
##'   \item{LoAge}{Lower age limit of the surveyed individuals in years.}
##'   \item{Ex}{Number of individuals examined for malaria.}
##'   \item{Pf}{Number of individuals tested positive for Plasmodium falciparum (malaria parasite).}
##'   \item{PfPR2.10}{Plasmodium falciparum parasite rate in the population (aged 2-10 years).}
##'   \item{Method}{Method used for malaria diagnosis (e.g., Rapid Diagnostic Test (RDT)).}
##'   \item{EVI}{Enhanced Vegetation Index, indicating vegetation density.}
##'   \item{Temperature}{Average temperature in degrees Celsius.}
##'   \item{Precipitation}{Total precipitation in millimeters.}
##'   \item{Population}{Population density in the area (in thousands).}
##'   \item{ITN}{Percentage of households with at least one insecticide-treated net (ITN).}
##'   \item{NTL}{Nighttime light intensity, indicating urbanization and infrastructure.}
##'   \item{Urban.Rural}{Indicator of whether the area is urban ('U') or rural ('R').}
##'   \item{elogit}{Empirical logit transformation of the positive cases - `elogit(Pf, Ex)`}
##'   \item{geometry}{Simple feature geometry (POINT) giving the sampling coordinates in metres}
##' }
##' @details The coordinate reference system is 32736.
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




