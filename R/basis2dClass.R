#'
setClass("basisfd")
#' A Class of 2-dimensional functional data objects
#'
#' The \code{basis2Dfd} class represents functional data ...
#' @slot bs1 a basisfd
#' @slot bs2  a basisfd ...
#' @slot support a matrix....
#'
#' @aliases basis2Dfd
#' @examples
#' x <- 1
setClass("basis2Dfd", slots = c(
  bs1 = "basisfd",
  bs2 = "basisfd",
  support = "matrix"
))

# Validity checks for basis2Dfd objects

setValidity("basis2Dfd", function(object) {
  bs1 <- object@bs1
  bs2 <- object@bs2
  return(TRUE)
})


#' Constructor for basis2Dfd objects
#'
#' @param bs1 a basisfd ...
#' @param bs2 a basisfd ...
#' 
#' @name basis2Dfd-constructor
#' @docType methods
#' @export  basis2Dfd
#' @keywords internal
#' @export
setGeneric("basis2Dfd", function(bs1,bs2){standardGeneric("basis2Dfd")})

#' @describeIn basis2Dfd Constructor for basis2Dfd objects when 
#' @param bs1 a basisfd ...
#' @param bs2 a basisfd ...
#' @docType methods
#' 
setMethod("basis2Dfd",
  signature = c(bs1 = "basisfd", bs2 = "basisfd"),
  function(bs1, bs2) {
    new("basis2Dfd", bs1 = bs1, bs2 = bs2, support = cbind(bs1$rangeval, bs2$rangeval))
  }
)
