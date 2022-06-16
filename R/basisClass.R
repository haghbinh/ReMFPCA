
#================================basismfd class=============================================#
# Class of multivariate functional data objects


#' A class for univariate functional data
#' 
#' The \code{basismfd} class represents functional data ...
#' @slot p a numeric ....
#' @slot support a list ....
#' @slot B a matrix ....
#' 
#' @aliases basismfd
#' 
#' @import methods
NULL

setClass("basismfd", slots = c(
  p = "numeric",
  support = "list",
  B = "matrix"
))

# Validity checks for basismfd objects
setValidity("basismfd", function(object) {
  p <- object@p
  s <- object@support
  if (length(p) > 1) {
    return("The argument p must have length 1.")
  }
  if (p != floor(p) | p < 0) {
    return("The argument p must be positive integer.")
  }
  return(TRUE)
})


#' Constructor for basismfd objects, third argument (B) passed as matrix or array of numerics
#' 
#' @param p a numeric specify number of variables
#' @param support a list ...
#' @param B a list ...
#' 
#' @name basismfd-constructor
#' @docType methods
#' @export basismfd
#' @keywords internal
#' 
setGeneric("basismfd", function(p, support, B){standardGeneric("basismfd")})



#' @describeIn basismfd Constructor for basismfd objects with \code{B} given as matrix.
#' @param p a numeric specify number of variables
#' @param support a list ...
#' @param B a matrix ...
#' @docType methods
#' 
setMethod("basismfd",
  signature = c(p = "numeric", support = "list", B = "matrix"),
  function(p, support, B) {
    new("basismfd", p = p, support = s, B = B)
  }
)