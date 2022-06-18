#================================basisEmp class=============================================#

#' A class for univariate Empirical functional data
#' 
#' @slot support a numeric ...
#' @slot grids a numeric ....
#' @slot B a matrix ....
#' 
#' @aliases basisEmp
#' 
#' @import methods
setClass("basisEmp", slots = c(
  support = "numeric",
  grids = "numeric",
  B = "matrix"
))

setValidity("basisEmp", function(object) {
  s <- object@support
  g <- object@grids
  B <- object@B
  if(length(g)!=nrow(B))
    return('Number of rows in B and length of grids must be equal.')
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
#' @param support a numeric ...
#' @param grids a numeric ...
#' @param B a matrix ...
#' @docType methods
#' 
setMethod("basisEmp",
          signature = c(support= "numeric", grids = "numeric", B = "matrix"),
          function(support, grids, B) {
            new("basisEmp",support= support, grids = grids, B = B)
          }
)





