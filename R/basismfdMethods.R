### ==============Support dimension method for the basismfd object =========
#' Support dimension of functional data
#'
#' This function returns the support dimension of an object of class  \code{basismfd}.
#' @param object An object of  class \code{mfd}.
#'  @return a numeric vector, giving the support dimension of each variable.
#' @export dimSupp
#' @examples
#' x <- 1
#' 
setGeneric("dimSupp", function(object) {standardGeneric("dimSupp")})

#' dimSupp for basismfd objects
#'
#' @keywords internal
setMethod("dimSupp", signature = "basismfd",
          function(object){
            s <- object@support
            vapply(s, FUN = ncol, FUN.VALUE = 0)})





### ==============Evaluation method for the basismfd object =========
#'  Values of basismfd Basis Functions 
#'
#' Computes the basis matrix evaluated at arguments in EVALARG associated
#' with basismfd object 

#' @param evalarg a list of numeric vectors of argument values at which the \code{basismfd} is to be evaluated.
#' @param basisobj  an object of class \code{basismfd} defining basis functions whose values are to be computed.
#'  @return a matrix of basis function values with rows corresponding to argument values and columns to basis functions.
#' @export eval.basismfd
#' @examples
#' x <- 1
#' 

setGeneric("eval.basismfd", function(evalarg, basisobj) {standardGeneric("eval.basismfd")})

#' evalvaluation of basismfd objects
#'
#' @keywords internal
setMethod("eval.basismfd",
  signature = c(evalarg = "list", basisobj = "basismfd"),
  function(evalarg, basisobj) {
    s <- basisobj@support
    p <- basisobj@p
    B <- basisobj@B
    d <- dimSupp(basisobj)
    if (length(evalarg) != p) {
      stop(paste("The length of evalarg must be ", p, " ."))
    }

    out <- list()
    for (i in 1:p) {
      out <- ifelse(is.basis(B[[i]]), eval.basis(evalarg[[i]]), B[[i]]@Basis)
    }
    return(out)
  }
)
