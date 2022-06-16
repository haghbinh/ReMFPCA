
#================================basismfd class=============================================#
#' Class of multivariate functional data objects

setClass("basismfd", slots = c(
  p = "numeric",
  support = "list",
  B = "list"
))

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



#' Constructor for basismfd objects
#' 
#' @param p a numeric specify number of variables
#' @param s a list ...
#' @param B a list ...
#' @name basismfd-constructor
#' @docType methods
#' @export basismfd
#' @keywords internal
#' 
#'  
setGeneric("basismfd", function(p, s, B) {
  standardGeneric("funData")
})


#' @describeIn basismfd Constructor for functional data objects with \code{argvals} given as list.
#' @param p a numeric specify number of variables
#' @param s a list ...
#' @param B a list ...
#' @docType methods
setMethod("basismfd",
  signature = c(p = "numeric", support = "list", B = "list"),
  function(argvals, X) {
    new("basismfd", p = p, support = s, B = B)
  }
)


#================================basisEmp class=============================================#


setClass("basisEmp", slots = c(
  grids = "numeric",
  B = "matrix"
))

setValidity("basisEmp", function(object) {
  g <- object@grids
  B <- object@B
  return(TRUE)
})

# Constructor of basisEmp class
Ebs <- function(g,B) {
  new("basisEmp", grids = g, B=B)
}



#================================basisfd class==============================================#


# Constructor of basisfd class
#' @importFrom  fda create.bspline.basis create.constant.basis create.exponential.basis create.fourier.basis create.monomial.basis create.polygonal.basis  create.power.basis
#' @export
create.basisfd <- function(rangeval = c(0, 1), nbasis = NULL, type = "bspline", ...) {
  out <- switch(type,
    "bspline" = create.bspline.basis(rangeval = rangeval, nbasis = nbasis, ...),
    "constant" = create.constant.basis(rangeval = rangeval, ...),
    "exponential" = create.exponential.basis(rangeval = rangeval, nbasis = nbasis, ...),
    "fourier" = create.fourier.basis(rangeval = rangeval, nbasis = nbasis, ...),
    "monomial" = create.monomial.basis(rangeval = rangeval, nbasis = nbasis, ...),
    "polygonal" = create.polygonal.basis(rangeval = rangeval, ...),
    "power" = create.power.basis(rangeval = rangeval, nbasis = nbasis, ...)
  )
}

# A = bs(1,2,c(2,5),list(2))


