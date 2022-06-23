#' A Class of Multidimensional Functional Data objects
#' @description The \code{mvfd} class represents functional data ...
#' @field basis a basismfd object
#' @field coeff a matrix with nrow=subjects and ncol=total number of basis ...
#' 
#' @examples
#' x <- 1
#' 
#' @importFrom fda is.basis eval.basis Data2fd
#' 
#' @export
mvfd <- R6::R6Class("mvfd",
  public = list(
    #' @description
    #' Constructor for mvfd objects
    #' @param mfd_list a list of mfd objects
    initialize = function(mfd_list) {
      init_mvfd_check(mfd_list)
      if (inherits(mfd_list, "mfd")) {
        mfd_list <- list(mfd_list)
      }
      private$.nobs <- mfd_list[[1]]$.nobs
      basis_list <- list()
      private$.nvar <- length(mfd_list)
      for (i in 1:private$.nvar) {
        mfd_list[[i]] <- mfd_list[[i]]$clone()
        basis_list[[i]] <- mfd_list$basis
        private$.coefs <- mfd_list[[i]]$coefs

        if (is.matrix(mfd_list[[i]]$coefs)) {
          private$.coefs <- rbind(private$.coefs, mfd_list[[i]]$coefs)
        } else {
          private$.coefs <- rbind(private$.coefs, apply(mfd_list[[i]]$coefs, 3, as.vector))
        }
      }
      private$.basis <- mvbasismfd$new(basis_list)
    },
    #' @description evalmfd
    #' @param evalarg a list of numeric vector of argument values at which the \code{mvfd} is to be evaluated.
    eval = function(evalarg) {
      eval_mvfd_validity_check(evalarg)
      eval_bs <- private$.basis$eval(evalarg)
      out <- list()
      for (i in 1:private$.nvar) {
        if (is.matrix(private$.coefs[[i]])) {
          out[[i]] <- eval_bs[[i]][[1]] %*% private$.coefs[[i]]
        } else {
          out[[i]] <- (eval_bs[[i]][[2]] %x% eval_bs[[i]][[1]]) %*% apply(private$.coefs[[i]], 3, as.vector)
        }
      }
      return(out)
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
    coefs = function(value) {
      if (missing(value)) {
        array(private$.coefs, c(unlist(mdbs2$nbasis), dim(X)[3]))
      } else {
        stop("`$coefs` is read only", call. = FALSE)
      }
    },
    nvar = function(value) {
      if (missing(value)) {
        private$.nvar
      } else {
        stop("`$nvar` is read only", call. = FALSE)
      }
    },
    nobs = function(value) {
      if (missing(value)) {
        private$.nobs
      } else {
        stop("`$nobs` is read only", call. = FALSE)
      }
    }
  ),
  private = list(
    .basis = NULL,
    .coefs = NULL, # we record vetorized of the coefs
    .nobs = NULL,
    .nvar = NULL
  )
)

# a function to check the validity of initializer
init_mvfd_check <- function(mfd_list) {
  if (inherits(mfd_list, "mfd")) {
    mfd_list <- list(mfd_list)
  }
  if (!all(sapply(mfd_list, function(x) inherits(x, "mfd")))) {
    stop("All the elements of the inputs list must have class of mfd")
  }
  n <- mfd_list[[1]]$nobs ;browser()
  for (y in mfd_list) {
    if (n != y$nobs) stop("Number of observations in all variables should be equal.")
  }
}
# a function to check the validity of evaluation
eval_mvfd_validity_check <- function(evalarg, dimSupp) {
  x <- 1
}

