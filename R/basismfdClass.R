#' @title A Class of Multidimensional Basis Functions
#' 
#' @description
#' The basismfd class represents functional data with multidimensional basis functions.
#' 
#'
#' @field basis A list of basis objects from the fda library.
#' @field dimSupp The dimension of the support domain of the `basismfd` object.
#' @field supp The matrix representing the ranges of the dimensions.
#' @field gram The Gram matrix.
#' @field nbasis A numeric vector containing the number of bases.
#' 
#' 
#' @importFrom fda is.basis eval.basis inprod
#' @importFrom Matrix Matrix
#' @export
basismfd <- R6::R6Class("basismfd",
  public = list(
    #' @description
        #'  Constructor for `basismfd` objects (same as Basismfd(...) )
    #' @usage Basismfd(...)
    #' @param ... A list of `basisfd` objects
    initialize = function(...) {
      basis <- list(...)
      init_basismfd_check(basis)
      
      if (is.basis(basis)) {
        private$.gram <- Matrix(inprod(basis, basis))
        private$.basis <- list(basis)
        private$.dimSupp <- 1
        private$.supp <- matrix(basis$rangeval, nrow = 2, ncol = 1)
        private$.nbasis <- basis$nbasis
      } else {
        private$.basis <- basis
        private$.dimSupp <- length(basis)
        private$.supp <- matrix(0, nrow = 2, ncol = length(basis))
        private$.gram <- 1
        for (i in 1:private$.dimSupp) {
          private$.gram <- Matrix(inprod(basis[[i]], basis[[i]]) %x% private$.gram)
          private$.nbasis[i] <- basis[[i]]$nbasis
          private$.supp[, i] <- basis[[i]]$rangeval
        }
      }
    },
    #' @description
      #' Evaluate the `basismfd` object at given argument values
    #' @param evalarg A list of numeric vectors of argument values at which the `basismfd` is to be evaluated
    #' @return A list of evaluated values
    eval = function(evalarg) {
      eval_basismf_validity_check(evalarg, private$.dimSupp)
      
      if (is.numeric(evalarg)) {
        evalarg <- list(evalarg)
      }
      
      out <- list()
      for (i in 1:length(evalarg)) {
        out[[i]] <- eval.basis(evalarg[[i]], private$.basis[[i]])
      }
      return(out)
    },
    #' @description
      #' Print method for `basismfd` objects
    #'
    #' @param ... Additional arguments to be passed to `print`
    print = function(...) {
      for (i in 1:private$.dimSupp) {
        cat("basis ", i, ":\n", sep = "")
        cat("type:", private$.basis[[i]]$type)
        cat("\nnbasis:", private$.nbasis[i])
        cat("\nsupport:", private$.supp[, i], "\n")
        if (i != private$.dimSupp) cat("\n")
      }
      invisible(self)
    }
  ),
  active = list(
    #' Getter and setter for `basis` field
    basis = function(value) {
      if (missing(value)) {
        private$.basis
      } else {
        stop("`$basis` is read only", call. = FALSE)
      }
    },
    
    #' Getter and setter for `dimSupp` field
    dimSupp = function(value) {
      if (missing(value)) {
        private$.dimSupp
      } else {
        stop("`$dimSupp` is read only", call. = FALSE)
      }
    },
    
    #' Getter and setter for `nbasis` field
    nbasis = function(value) {
      if (missing(value)) {
        private$.nbasis
      } else {
        stop("`$nbasis` is read only", call. = FALSE)
      }
    },
    
    #' Getter and setter for `supp` field
    supp = function(value) {
      if (missing(value)) {
        private$.supp
      } else {
        stop("`$supp` is read only", call. = FALSE)
      }
    },
    
    #' Getter and setter for `gram` field
    gram = function(value) {
      if (missing(value)) {
        private$.gram
      } else {
        stop("`$gram` is read only", call. = FALSE)
      }
    }
  ),
  
  private = list(
    .basis = NULL,
    .dimSupp = NULL,
    .nbasis = NULL,
    .supp = NULL,
    .gram = NULL
  )
)

#' @export
Basismfd <- function(...) {
  basismfd$new(...)
}

#' Check the validity of the initializer
#' @param basis A list of `basisfd` objects
init_basismfd_check <- function(basis) {
  if (!is.basis(basis) & is.list(basis)) {
    if (!all(sapply(basis, is.basis))) {
      stop("All elements of the basis list must be `basisfd` objects.")
    }
  }
}

#' Check the validity of the evaluation
#' @param evalarg A list or numeric vector
#' @param dimSupp Dimension of the `basismfd` object
eval_basismf_validity_check <- function(evalarg, dimSupp) {
  if (!is.list(evalarg) & !is.numeric(evalarg)) {
    stop("evalarg must be a list or numeric vector.")
  }
  if (!all(sapply(evalarg, is.numeric))) {
    stop("evalarg must be a list of numeric vectors.")
  }
  if (is.numeric(evalarg)) {
    evalarg <- list(evalarg)
  }
  if (length(evalarg) != dimSupp) {
    stop("Length of evalarg list must be equal to dimSupp.")
  }
}
