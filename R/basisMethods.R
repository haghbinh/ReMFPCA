### ==============Support dimension method for the basismfd object =========
#' Support dimension of functional data
#' 
#' This function returns the support dimension of an object of class  \code{basismfd}.
#' @param object An object of  class \code{mfd}.
#'  @return a numeric vector, giving the support dimension of each variable.
#' @export
setGeneric("dimSupp", function(object) {standardGeneric("dimSupp")})
setMethod("dimSupp", signature = "basismfd",
          function(object){
            s <- object@support
            vapply(s, FUN = ncol, FUN.VALUE = 0)})


### ==============Eval method for the basismfd object =========
