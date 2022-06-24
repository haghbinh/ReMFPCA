#' A Class of Multidimensional Functional Data objects
#' @description The \code{mfd} class represents functional data ...
#' @field basis a basismfd object
#' @field coeff a matrix with nrow=subjects and ncol=total number of basis ...
#'
#' @examples
#' x <- 1
#' @importFrom fda is.basis eval.basis Data2fd
#'
#' @export
mfd <- R6::R6Class("mfd",
  public = list(
    #' @description
    #' Constructor for mfd objects
    #' @param mdbs a basismfd object
    initialize = function(argval = NULL, X, mdbs, method = c("data", "coefs")) {
      init_mfd_check(argval, X, mdbs, method)
      if (is.basis(mdbs)) {
        mdbs <- basismfd$new(mdbs)
      } else {
        mdbs <- mdbs$clone()
      }
      private$.basis <- mdbs
      if (is.vector(X)) X <- matrix(X)
      if (method[1] == "coefs") {
        if (length(dim(X))>2) X <- apply(X, length(dim(X)), as.vector)
        else if (mdbs$dimSupp>1 && ncol(X)!=1) X <- matrix(X)
        private$.coefs <- X
      } else {
        if (is.null(argval) && method != "coefs") {
          argval <- list()
          for (j in 1:mdbs$dimSupp) {
            argval[[j]] <- seq(mdbs$supp[1, j], mdbs$supp[2, j], len = dim(X)[j])
          }
        }
        if (is.numeric(argval)) argval <- list(argval)

        Bmat <- private$.basis$eval(argval)
        if (length(argval) > 1) { # This is for 2D case
          B <- Bmat[[2]] %x% Bmat[[1]]
          if (is.matrix(X)) X <- array(X, dim = c(dim(X), 1))
          private$.coefs <- solve(t(B) %*% B) %*% t(B) %*% apply(X, 3, as.vector)
        } else { # This is for 1-dimensional case
          B <- Bmat[[1]]
          private$.coefs <- solve(t(B) %*% B) %*% t(B) %*% X
        }
      }
      private$.nobs <- tail(dim(X), 1)
    },
    #' @description evalmfd
    #' @param evalarg a list of numeric vector of argument values at which the \code{mfd} is to be evaluated.
    eval = function(evalarg) {
      eval_mfd_validity_check(evalarg)
      if (is.numeric(evalarg)) evalarg <- list(evalarg)
      Bmat <- private$.basis$eval(evalarg)
      if (length(evalarg) > 1) {
        Xhat <- Bmat[[2]] %x% Bmat[[1]] %*% private$.coefs
        Xhat <- array(Xhat, dim = c(sapply(evalarg, length), private$.nobs))
      } else {
        Xhat <- Bmat[[1]] %*% private$.coefs
      }
      return(Xhat)
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
        array(private$.coefs, c(unlist(private$.basis$nbasis), private$.nobs))
      } else {
        stop("`$coefs` is read only", call. = FALSE)
      }
    },
    nobs = function(value) {
      if (missing(value)) {
        private$.nobs
      } else {
        stop("`$coefs` is read only", call. = FALSE)
      }
    }
  ),
  private = list(
    .basis = NULL,
    .coefs = NULL, # we record vetorized of the coefs
    .nobs = NULL
  )
)

# a function to check the validity of initializer
init_mfd_check <- function(argval, X, basis, method) {
  x <- 1
}
# a function to check the validity of evaluation
eval_mfd_validity_check <- function(evalarg, dimSupp) {
  x <- 1
}
