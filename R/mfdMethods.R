### generic functions ###


### ==============Support dimension method for the mfd object (by Hossein)=========
#' Support dimension of functional data
#' 
#' This function returns the support dimension of an object of class  \code{mfd}.
#' @param object An object of  class \code{mfd}.
#'  @return a numeric vector, giving the support dimension of each variable.
#' @export dimSupp
setGeneric("dimSupp", function(object) {standardGeneric("dimSupp")})
setMethod("dimSupp", signature = "mfd",
          function(object){
            vapply(object@grid, FUN = function(x) ncol(as.matrix(x)), FUN.VALUE = 0)})
            

### ============== number of variables method for the mfd object (by Hossein)=========
#' Get the number of variables
#' 
#' This functions returns the number of variables in a \code{mfd} object.
#' @param object An object of class \code{mfd}
#' @return The number of variables in \code{mfd}.
#' @export nVar
setGeneric("nVar", function(object) {standardGeneric("nVar")})
setMethod("nVar", signature = "mfd",
          function(object){length(object@C)})


### ============== number of observations method for the mfd object (by Hossein)=========

#' Get the number of observations
#' 
#' This functions returns the number of observations in a \code{mfd} object.
#' @param object An object of class \code{mfd}
#' @return The number of observations in \code{mfd}.
#' @export nObs
setGeneric("nObs", function(object) {standardGeneric("nObs")})
setMethod("nObs", signature = "mfd",
          function(object){ncol(object@C[[1]])})


### ==============Show method for the mfd object (by Hossein)======================
#' A print method for univariate functional data
#' 
#' This function prints basic information about a \code{mfd} object. This is
#' the standard console output for \code{mfd} objects.
#' @export
#' 
print.mfd <- function(object,...){
  N <- nObs(object)
  p <- nVar(object)
  if(p==1){
    cat("Functional data with ", N ," observations of ", dimSupp(object) ,"-dimensional support.\n", sep = "")
    print('Grids:')
    print(object@grid)
    print("Coefficients:")
    print(object@C)
    print("Basis evaluated matrix:")
    print(object@B)
  } else {
    cat("Multivariate functional data with ", N ," observations from the ", p  ," variable vectors of (", paste(dimSupp(object),collapse = ',') ,")-dimensional support.\n", sep = "")
    for(i in 1:p){
      cat('Grids of the variable ',i,':\n')
      print(object@grid)
      cat("Coefficients of the variable ",i,':\n')
      print(object@C)
      cat("Basis evaluated matrix of the variable ",i,':\n')
      print(object@B)
    }
  }
}


setMethod("show", signature = "mfd", function(object) {print.mfd(object)})




### ============== evaluation of mfd object (by Hossein)=========


### =====================Plot methods for mft objects================== ###
#' @export
plot.mfd <- function(obj, vars = NULL, types = NULL, subplot = TRUE, mains = NULL, ylabels = NULL, xlabels = NULL, tlabels = NULL, zlabels = NULL, ...) {
  p <- nVar()
  N <- ncol(x@C[[1]])
}


