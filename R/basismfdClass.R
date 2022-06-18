#' @include basisEmpClass.R
#' @include basis2dClass.R
setClass("basisfd")
#================================basismfd class=============================================#
#' A Class of multivariate functional data objects
#'
#' The \code{basismfd} class represents functional data ...
#' @slot p a numeric....
#' @slot support a matrix ...
#' @slot B a list ....
#'
#' @aliases basismfd
#' @examples
#' require(fda)
#' require(ReMFPCA)
#' bs_fd <- create.fourier.basis(c(0,2*pi),5)
#' argval <- seq(0,2*pi,1)
#' b <- eval.basis(argval,bs_fd)
#' bs_Em <- basisEmp(c(0,2*pi),argval,b)
#' B <- list(bs_fd,bs_Em)
#' 
#' bs1 <- basismfd(B=bs_Em)
#' bs1@support
#' 
#' bs2 <- basismfd(bs_fd)
#' bs2@support
#' bs3 <- basismfd(B)
#' bs3@support
#' @importFrom fda is.basis 
setClass("basismfd", slots = c(
  p = "numeric",
  support = "list",
  B = "list"
))

# Validity checks for basismfd objects

setValidity("basismfd", function(object) {
  p <- object@p
  B <- object@B
  if (length(p)>1)
    return('p must be a numeric of length one.')
  if(is.list(B) & !is.basis(B) & !all(sapply(B, function(x) is.basis(x) | inherits(x,"basisEmp")| inherits(x,"basis2Dfd"))))
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
#' @param B a list ...
#' @docType methods
#' @export 
setMethod("basismfd",
  signature = c(B = "list"),
  function(B) {
    p <- length(B)
    new("basismfd",p = length(B),
      support = lapply(B, function(x) if (is.basis(x)) {
        return(as.matrix(x$rangeval))
      } else {
        return(as.matrix(x@support))
      }),
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
            new("basismfd",p = 1 ,support =list(B@support), B = list(B))
          }
)


#' @describeIn basismfd Constructor for basismfd objects nd \code{B} is given as basisfd(from fda package)
#'
#' @docType methods
setMethod("basismfd",
          signature = c( B = "basisfd"),
          function(B) {
            new("basismfd",p = 1 ,support =list(as.matrix(B$rangeval)), B = list(B))
          }
)



#' @describeIn basismfd Constructor for basismfd objects
#'
#' @docType methods
setMethod("basismfd",
          signature = c( B = "basis2Dfd"),
          function(B) {
            new("basismfd",p = 1 ,support =list(B@support), B = list(B))
          }
)

