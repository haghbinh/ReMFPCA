### ==============Evaluation method for the basismfd object =========
#'  Values of basismfd Basis Functions 
#'
#' Computes the basis matrix evaluated at arguments in EVALARG associated with basismfd object 
#' @param evalarg a list of numeric vectors of argument values at which the \code{basismfd} is to be evaluated.
#' @param basisobj  an object of class \code{basismfd} defining basis functions whose values are to be computed.
#' @return a matrix of basis function values with rows corresponding to argument values and columns to basis functions.
#' @export evalbasismfd
#' @examples
#' x <- 1
#' 
setGeneric("evalbasismfd", function(evalarg, basisobj) {standardGeneric("evalbasismfd")})

#' evalbasismfd of basismfd objects
#' @importFrom fda eval.basis
#' @keywords internal
setMethod("evalbasismfd",
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
              if(is.basis(B[[i]])){
               out[[i]] <- eval.basis(evalarg[[i]],B[[i]])
              }else{
                if(inherits(B[[i]],"basisEmp")){
                  out[[i]] <- B[[i]]@Basis
                } else{
                  b1 <- eval.basis(evalarg[[i]][[1]],B[[i]]@bs1)
                  b2 <- eval.basis(evalarg[[i]][[2]],B[[i]]@bs2)
                  out[[i]] <- kronecker(b1,b2)
                }
                
              }
            }
            return(out)
          }
)




#' evalbasismfd of basismfd objects
#' @importFrom fda eval.basis
#' @keywords internal
setMethod("evalbasismfd",
          signature = c(evalarg = "numeric", basisobj = "basismfd"),
          function(evalarg, basisobj) {
            s <- basisobj@support
            p <- basisobj@p
            B <- basisobj@B
            d <- dimSupp(basisobj)
            if (p > 1) {
              warning("The argument vector is applid for all basises")
            }
            out <- list()
            for (i in 1:p) {
              if(is.basis(B[[i]])){
                out[[i]] <- eval.basis(evalarg,B[[i]])
              }else{
                if(inherits(B[[i]],"basisEmp")){
                  out[[i]] <- B[[i]]@Basis
                } else{
                  b1 <- eval.basis(evalarg,B[[i]]@bs1)
                  b2 <- eval.basis(evalarg,B[[i]]@bs2)
                  out[[i]] <- kronecker(b1,b2)
                }
                
              }
            }
            return(out)
          }
)


#' Support dimension of functional data
#'
#' This function returns the support dimension of an object of class  \code{basismfd}.
#' @param object An object of  class \code{mfd}.
#' @return a numeric vector, giving the support dimension of each variable.
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


