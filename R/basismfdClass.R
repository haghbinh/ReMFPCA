
#================================basismfd class=============================================#
#' A Class of multivariate functional data objects
#' 
#' The \code{basismfd} class represents functional data ...
#' @slot p a numeric....
#' @slot B a list ....
#' 
#' @aliases basismfd
#' 
#' @import methods


setClass("basismfd", slots = c(
  p = "numeric",
  B = "list"
))

# Validity checks for basismfd objects
#' @importFrom fda is.basis
setValidity("basismfd", function(object) {
  p <- object@p
  B <- object@B
  if (length(p)>1)
    return('p must be a numeric of length one.')
  if(is.list(B) & !is.basis(B) & !all(sapply(B, function(x) is.basis(x) | inherits(x,"basisEmp"))))
    return('All the elements of B must have the classes of types basisfd or basisEmp.')
  return(TRUE)
})


#' Constructor for basismfd objects
#' 
#' @param B a list of basisEmp or basisfd 
#' 
#' @name basismfd-constructor
#' @docType methods
#' @export basismfd
#' @keywords internal
#' 
setGeneric("basismfd", function(B){standardGeneric("basismfd")})

#' @describeIn basismfd Constructor for basismfd objects when both \code{support} and \code{B} are given as list
#' @param support a list ...
#' @param B a list ...
#' @docType methods
#' 
setMethod("basismfd",
  signature = c( B = "list"),
  function(B) {
    new("basismfd", p = length(B), B = B)
  }
)


#' @describeIn basismfd Constructor for basismfd objects when \code{support} is given as vector
#' and \code{B} is given as matrix.
#'   
#' @docType methods
setMethod("basismfd",
          signature = c(B = "basisEmp"),
          function(B) {
            new("basismfd", p = 1, B = list(B))
          }
)


#' @describeIn basismfd Constructor for basismfd objects nd \code{B} is given as basisfd(from fda package)
#'   
#' @docType methods
setMethod("basismfd",
          signature = c( B = "basisfd"),
          function(B) {
            new("basismfd", p = 1, B = list(B))
          }
)



