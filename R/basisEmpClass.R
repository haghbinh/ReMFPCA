#================================basisEmp class=============================================#

#' A class for univariate Empirical functional data
#' 
#' @slot support a numeric ...
#' @slot grids a numeric ....
#' @slot Basis a matrix ....
#' 
#' @aliases basisEmp
#' 
#' @import methods
setClass("basisEmp", slots = c(
  support = "matrix",
  grids = "list",
  Basis = "matrix"
))

setValidity("basisEmp", function(object) {
  s <- object@support
  g <- object@grids
  B <- object@Basis
  # if(length(g)!=nrow(B))
  #   return('Number of rows in B and length of grids must be equal.')
  return(TRUE)
})




#' Constructor for basisEmp objects,
#' @param support a numeric ...
#' @param grids a numeric vector  ...
#' @param Basis a numeric matrix ... 
#' 
#' @name basisEmp-constructor
#' @docType methods
#' @export basisEmp
#' @keywords internal
#' 
setGeneric("basisEmp", function(support, grids, Basis){standardGeneric("basisEmp")})

#' @describeIn basisEmp Constructor for basisEmp objects when .....
#' @param support a numeric ...
#' @param grids a numeric ...
#' @param Basis a matrix ...
#' @docType methods
#' 
setMethod("basisEmp",
          signature = c(support= "matrix", grids = "list", Basis = "matrix"),
          function(support, grids, Basis) {
            new("basisEmp",support= support, grids = grids, Basis = Basis)
          }
)


#' @describeIn basisEmp Constructor for basisEmp objects .....
#'
#' @docType methods
setMethod("basisEmp",
          signature = c(support= "numeric", grids = "numeric", Basis = "matrix"),
          function(support, grids, Basis) {
            new("basisEmp",support= as.matrix(support), grids = list(grids), Basis = Basis)
          }
)

