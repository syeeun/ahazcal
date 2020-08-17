#' @title Sample Data for 'ahazcal' Package
#' @description Sample cohort with simulated nested case-control study where N = 1000 and 173 cases with m = 1 control per case. This dataset is used to run the provided examples for the functions in the package 'ahazcal'.
#' @format A data frame with N = 1000 rows and 10 variables:
#' \describe{
#'   \item{eventime}{time-to-event (event = either case or censored)}
#'   \item{ind.fail}{binary outcome status (1 = case; 0 = non-case)}
#'   \item{X1}{binary model covariates with time-varying effects; available in a full cohort}
#'   \item{X2}{binary model covariates with time-varying effects; available in a full cohort}
#'   \item{Z1}{binary model covariates with time-invariant effects; available in a full cohort}
#'   \item{Z2}{binary model covariates with time-invariant effects; only available in nested case-control samples; NA if missing}
#'   \item{U}{an ancillary predictor for Z2} 
#'   \item{ind.ph2}{phase 2 inclusion status (1 = included as a case or selected control; 0 = otherwise)}
#'   \item{nrisk}{number of at-risk subjects at each 'eventime'}
#'   \item{incl.prob}{inclusion probability }
#'}
"sample_data"