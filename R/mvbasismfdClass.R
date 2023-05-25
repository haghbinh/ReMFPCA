#' @title A Class of Multivariate Multidimensional Basis Functions
#'
#' @description
#' The `mvbasismfd` class represents functional data ...
#'
#' An object of class "mvbasismfd" is a list containing the following elements:
#' @field nvar number of variables
#' @field basis A list of mvbasisfd objects
#' @field dimSupp A sequence of positive integers specifying the dimension...
#' @field nbasis A list of integers specifying the number of basis functions...
#' @field supp A list of matrices specifying the support of basis functions...
#' @field gram A block diagonal matrix
#'
#
#'
#' @importFrom fda is.basis eval.basis
#' @importFrom Matrix Matrix bdiag
#' @seealso \code{\link{mvmfd}}, \code{\link{basismfd}}
#' @export
mvbasismfd <- R6::R6Class("mvbasismfd",
  public = list(
    #' @description
    #' Constructor for `mvbasismfd` objects (same as Mvbasismfd(...) )
    #' 
    #' @usage Mvbasismfd(basis)
    #' 
    #' @param basis A list of `basismfd` objects
    initialize = function(basis) {
      if (is.basis(basis) | is.basismfd(basis)) basis <- list(basis)
      init_mvbasismfd_check(basis)
      private$.nvar <- length(basis)
      private$.gram <- bdiag()
      for (i in 1:private$.nvar) {
        if (is.basis(basis[[i]])) {
          basis[[i]] <- basismfd$new(basis[[i]])
        } else {
          basis[[i]] <- basis[[i]]$clone()
        }
        private$.dimSupp[i] <- basis[[i]]$dimSupp
        private$.nbasis[[i]] <- basis[[i]]$nbasis
        private$.supp[[i]] <- basis[[i]]$supp
        private$.gram <- bdiag(private$.gram, basis[[i]]$gram)
      }
      private$.basis <- basis
    },

    #' @description
    #' Evaluate the `mvbasismfd` object at given argument values
    #'
    #' @param evalarg A list of numeric vectors of argument values at which the `basismfd` is to be evaluated
    #' @return A list of evaluated values
    eval = function(evalarg) {
      eval_mvbasismf_validity_check(evalarg, private$.nvar)
      if (is.numeric(evalarg)) {
        evalarg <- list(list(evalarg))
      }
      out <- list()
      for (i in 1:length(evalarg)) {
        out[[i]] <- private$.basis[[i]]$eval(evalarg[[i]])
      }
      return(out)
    }
  ),
  active = list(
    # nvar field
    nvar = function(value) {
      if (missing(value)) {
        private$.nvar
      } else {
        stop("`$nvar` is read only", call. = FALSE)
      }
    },

    # basis field
    basis = function(value) {
      if (missing(value)) {
        private$.basis
      } else {
        stop("`$basis` is read only", call. = FALSE)
      }
    },

    # dimSupp field
    dimSupp = function(value) {
      if (missing(value)) {
        private$.dimSupp
      } else {
        stop("`$dimSupp` is read only", call. = FALSE)
      }
    },

    # nbasis field
    nbasis = function(value) {
      if (missing(value)) {
        private$.nbasis
      } else {
        stop("`$nbasis` is read only", call. = FALSE)
      }
    },

    # supp field
    supp = function(value) {
      if (missing(value)) {
        private$.supp
      } else {
        stop("`$supp` is read only", call. = FALSE)
      }
    },

    # gram field
    gram = function(value) {
      if (missing(value)) {
        private$.gram
      } else {
        stop("`$gram` is read only", call. = FALSE)
      }
    }
  ),
  private = list(
    .nvar = NULL,
    .basis = NULL,
    .dimSupp = NULL,
    .nbasis = list(),
    .supp = list(),
    .gram = NULL
  )
)
#' @rdname mvbasismfd
#' @seealso \code{\link{basismfd}}

#' @title Constructor for mvbasismfd objects
#'
#' @description
#' Constructor for `mvbasismfd` objects (same as Mvbasismfd(...) )
#'
#' @import R6
#' @param basis A list of basisfd objects
#' @export
Mvbasismfd <- function(basis) mvbasismfd$new(basis)
#' @rdname mvbasismfd
#' @seealso \code{\link{basismfd}}

#' 
#' @param mvbasismfd_obj An 'mvmfd' object
#' @param i An index or indices specifying the subsets to extract for the first dimension
#' @return An 'mvbasismfd' object containing the specified subsets
#' @export
"[.mvbasismfd" <- function(mvbasismfd_obj, i = "index") {
  if(max(i)> mvbasismfd_obj$nvar | min(i)<1) stop(" subscript i out of bounds")
  if(length(i)==1){
    return(mvbasismfd_obj$basis[[i]])
  }else{
    mvbasismfd_list <- list()
    for(j in 1:length(i)) {
      mvbasismfd_list[[j]] <- mvbasismfd_obj$basis[[j]]
    }    
    return(mvbasismfd$new(mvbasismfd_list))
  }
}
