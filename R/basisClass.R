
#================================basismfd class=============================================#
#' A Class of multivariate functional data objects
#' 
#' The \code{basismfd} class represents functional data ...
#' @slot p a numeric....
#' @slot support a matrix ...
#' @slot B a list ....
#' 
#' @aliases basismfd
#' 
#' @import methods


setClass("basismfd", slots = c(
  p = "numeric",
  support = "list",
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
  signature = c(B = "list"),
  function(B) {
    p <- length(B)
    new("basismfd",
      support = lapply(B, function(x) if (is.basis(x)) {
        return(as.matrix(x$rangeval))
      } else {
        return(as.matrix(x@support))
      }),
      p = length(B),
      B = B
    )
  }
)


#' @describeIn basismfd Constructor for basismfd objects when \code{support} is given as vector
#' and \code{B} is given as matrix.
#'   
#' @docType methods
setMethod("basismfd",
          signature = c(B = "basisEmp"),
          function(B) {
            new("basismfd",support =list(B@support),  p = 1, B = list(B))
          }
)


#' @describeIn basismfd Constructor for basismfd objects nd \code{B} is given as basisfd(from fda package)
#'   
#' @docType methods
setMethod("basismfd",
          signature = c( B = "basisfd"),
          function(B) {
            new("basismfd",support =list(as.matrix(B$rangeval)) , p = 1, B = list(B))
          }
)






#================================basisEmp class=============================================#

#' A class for univariate Empirical functional data
#' 
#' @slot support a matrix ...
#' @slot grids a numeric ....
#' @slot B a matrix ....
#' 
#' @aliases basisEmp
#' 
#' @import methods
setClass("basisEmp", slots = c(
  support = "matrix",
  grids = "list",
  B = "matrix"
))

setValidity("basisEmp", function(object) {
  s <- object@support
  g <- object@grids
  B <- object@B
  # if(length(g)!=nrow(B))
  #   return('Number of rows in B and length of grids must be equal.')
  return(TRUE)
})



#' Constructor for basisEmp objects,
#' @param support a numeric ...
#' @param grids a numeric vector  ...
#' @param B a numeric matrix ... 
#' 
#' @name basisEmp-constructor
#' @docType methods
#' @export basisEmp
#' @keywords internal
#' 
setGeneric("basisEmp", function(support, grids, B){standardGeneric("basisEmp")})

#' @describeIn basisEmp Constructor for basisEmp objects when .....
#' @param support a 2*k matrix  where the first row specify min and second max and k is dimension...
#' @param grids a list of length k  ...
#' @param B a matrix ...
#' @docType methods
#' 
setMethod("basisEmp",
          signature = c(support= "matrix", grids = "list", B = "matrix"),
          function(support, grids, B) {
            new("basisEmp",support= support, grids = grids, B = B)
          }
)



setMethod("basisEmp",
          signature = c(support= "numeric", grids = "numeric", B = "matrix"),
          function(support, grids, B) {
            new("basisEmp",support= as.matrix(support), grids = list(grids), B = B)
          }
)

