#' Electrical Power Dataset
#'
#' This dataset contains measurements of voltage and electrical power consumption
#' recorded from a household between December 2006 and November 2010. 
#'
#' @name electrical_power_data
#' @docType data
#' @format A bivariate functional data object of class `mvmfd` with the following fields:
#' \describe{
#'   \item{voltage}{Voltage measurements in volts.}
#'   \item{power}{Calculated power consumption in watts.}
#' }
#' @source The dataset was collected from https://weather.com.
#' @examples
#' \dontrun{
#' # Load the Electrical Power Dataset
#' data("electrical_power_data")
#' head(electrical_power_data)
#' }
NULL
