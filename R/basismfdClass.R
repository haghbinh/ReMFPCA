#' A Class of Multidimensional Basis Functions
#' @description The \code{basismfd} class represents functional data ...
#' @field basis a list of basisfd objects
#' @field dimSupp  a positive integer specify the dimension...
#' 
#' @examples
#' x <- 1
#' 
#' @importFrom fda is.basis eval.basis inprod
#' @importFrom Matrix Matrix
#' 
#' @export
#' 
basismfd <- R6::R6Class("basismfd",
  public = list(
    #' @description
    #' Constructor for basismfd objects
    #' @param basis a list of basisfd objects
    initialize = function(basis) {
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
    #' @description evalbasismfd
    #' @param evalarg a list of numeric vector of argument values at which the \code{basismfd} is to be evaluated.
    #' @return a list
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
    print = function(...) {
      for(i in 1:private$.dimSupp){
        cat("basis ",i,":\n",sep = "")
        cat("type:",private$.basis[[i]]$type)
        cat("\nnbasis:",private$.nbasis[i])
        cat("\nsupport:",private$.supp[,i],"\n")
        if (i!=private$.dimSupp) cat("\n")
      }
      invisible(self)
    }
  ),
  active = list(
    basis = function(value) {
      if (missing(value)) {
        private$.basis
      } else {
        stop("`$basis` is read only", call. = FALSE)
      }
    },
    dimSupp = function(value) {
      if (missing(value)) {
        private$.dimSupp
      } else {
        stop("`$dimSupp` is read only", call. = FALSE)
      }
    },
    nbasis = function(value) {
      if (missing(value)) {
        private$.nbasis
      } else {
        stop("`$nbasis` is read only", call. = FALSE)
      }
    },
    supp = function(value) {
      if (missing(value)) {
        private$.supp
      } else {
        stop("`$supp` is read only", call. = FALSE)
      }
    },
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

# a function to check the validity of initializer
init_basismfd_check <- function(basis) {
  if (!is.basis(basis) & is.list(basis)) {
    if (!all(sapply(basis, is.basis))) {
      stop("All the elements of basis list must be basisfd object.")
    }
  }
}
# a function to check the validity of evaluation
eval_basismf_validity_check <- function(evalarg, dimSupp) {
  if (!is.list(evalarg) & !is.numeric(evalarg)) {
    stop("evalarg must be a list or numeric vector")
  }
  if (!all(sapply(evalarg, is.numeric))) {
    stop("evalarg must be a list of numeric vector")
  }
  if (is.numeric(evalarg)) {
    evalarg <- list(evalarg)
  }
  if (length(evalarg) != dimSupp) {
    stop(" length of evalarg list must be equal to dimSupp")
  }
}

#' @export
Basismfd <- function(basis) basismfd$new(basis)
