
#================================basismfd class=============================================#
# Class of multivariate functional data objects


#' A class for univariate functional data
#' 
#' The \code{basismfd} class represents functional data ...
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
  B = "list"
))

# Validity checks for basismfd objects
#' @importFrom fda is.basis
setValidity("basismfd", function(object) {
  p <- object@p
  s <- object@support
  B <- object@B
  if (length(p)>1)
    return('p must be a numeric of length one.')
  if(is.list(B) & !is.basis(B) & is.numeric(s))
    return('support must be a list')
  if(is.numeric(s) & length(s) != 2)
    return('When the support is vector, the length must be 2.')
  return(TRUE)
  
})


#' Constructor for basismfd objects, third argument (B) passed as matrix or array of numerics
#' 
#' @param support a numeric vector or list ...
#' @param B a numeric matrix or list 
#' 
#' @name basismfd-constructor
#' @docType methods
#' @export basismfd
#' @keywords internal
#' 
setGeneric("basismfd", function(support, B){standardGeneric("basismfd")})



#' @describeIn basismfd Constructor for basismfd objects when both \code{support} and \code{B} are given as list
#' @param support a list ...
#' @param B a list ...
#' @docType methods
#' 
setMethod("basismfd",
  signature = c(support = "list", B = "list"),
  function(support, B) {
    new("basismfd", p = length(B), support = support, B = B)
  }
)


#' @describeIn basismfd Constructor for basismfd objects when \code{support} is given as vector
#' and \code{B} is given as matrix.
#'   
#' @docType methods
setMethod("basismfd",
          signature = c(support = "numeric", B = "matrix"),
          function(support, B) {
            new("basismfd", p = 1, support = list(support), B = list(B))
          }
)


#' @describeIn basismfd Constructor for basismfd objects nd \code{B} is given as basisfd(from fda package)
#'   
#' @docType methods
setMethod("basismfd",
          signature = c(B = "basisfd"),
          function(support, B) {
            new("basismfd", p = 1, support = list(B$rangeval), B = list(B))
          }
)



