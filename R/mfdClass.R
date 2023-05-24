#' @title  A Class of Multidimensional Functional Data objects
#'
#' @description
#' The `mfd` class represents class of multidimensional functional data.
#'
#' @field basis an object of the class basismfd or fda::bsaia,
#' @field coefs  record vetorized of the coefs
#' @field nobs number of the observation
#'
#' @examples
#' require(fda)
#' bs1 <- create.fourier.basis(c(0,2*pi),5)
#' bs2 <- create.bspline.basis(c(0,1),7)
#' bs3 <- create.exponential.basis(c(0,2),3)
#' 
#' #1-D mfd :_____________________________
#' argval <- seq(0,2*pi,length.out=100)
#' nobs <- 10;
#' X <- outer(sin(argval),seq(0.5,1.5,length.out=nobs))
#' mdbs1 <- Basismfd(bs1)
#' mfd1 <- Mfd(X=X, mdbs = mdbs1)
#' inprod_mfd(mfd1,mfd1)
#' norm_mfd(mfd1)
#' mfd0 <- 2.5*mfd1
#' mfd1-mfd0
#' mfd1[1:3]
#' 

#' mfd1$eval(argval)
#' mfd1c <- Mfd(X=mfd1$coefs, mdbs = mdbs1, method = "coefs")
#' all.equal(c(mfd1$basis,mfd1$coefs,mfd1$nobs),c(mfd1c$basis,mfd1c$coefs,mfd1c$nobs))
#' length(mfd1)
#' mean(mfd1)
#' plot(mfd1)
#' 
#'
#' @importFrom fda is.basis eval.basis Data2fd
#'
#' @export
mfd <- R6::R6Class("mfd",
  public = list(
    #' @param argval A list of numeric vectors of argument values at which the `mfd` object is to be evaluated
    #' @param X A numeric matrix corresponds to basis expansion coefficients
    #' if `method="coefs"` and discrete observations if `method="data"`.
    #' @param mdbs a basismfd object
    #' @param method determine the `X` matrix type as "coefs" and "data".
    initialize = function(argval = NULL, X, mdbs, method = "data") { # c("data", "coefs")
      init_mfd_check(argval, X, mdbs, method)
      if (is.basis(mdbs)) {
        mdbs <- basismfd$new(mdbs)
      } else {
        mdbs <- mdbs$clone()
      }
      private$.basis <- mdbs
      if (is.vector(X)) X <- matrix(X)
      if (method[1] == "coefs") {
        if (length(dim(X)) > 2) {
          X <- apply(X, length(dim(X)), as.vector)
        } else if (mdbs$dimSupp > 1 && ncol(X) != 1) X <- matrix(X)
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
    #' @description Evaluation an `mfd` object in some arguments.
    #' @param evalarg a list of numeric vector of argument values at which the \code{mfd} is to be evaluated.
    #' @return A matrix of evaluated values
    eval = function(evalarg) {
      if (is.numeric(evalarg)) evalarg <- list(evalarg)
      Bmat <- private$.basis$eval(evalarg)
      if (length(evalarg) > 1) {
        Xhat <- Bmat[[2]] %x% Bmat[[1]] %*% private$.coefs
        Xhat <- array(Xhat, dim = c(sapply(evalarg, length), private$.nobs))
      } else {
        Xhat <- Bmat[[1]] %*% private$.coefs
      }
      return(Xhat)
    },
    #' @description
        #'  Print method for `mfd` objects
    #'
    #' @param ... Additional arguments to be passed to `print`
    print = function(...) {
      cat("A ", private$.basis$dimSupp, "-Dimensional 'mfd' object:", sep = "")
      cat("\nnobs:", private$.nobs, "\n")
      print(private$.basis)

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
    #' Getter and setter for `coefs` field
    coefs = function(value) {
      if (missing(value)) {
        array(private$.coefs, c(unlist(private$.basis$nbasis), private$.nobs))
      } else {
        stop("`$coefs` is read only", call. = FALSE)
      }
    },
    #' Getter and setter for `nobs` field
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
#' @rdname mfd
#' @seealso \code{\link{basismfd}}

#' @title  A Class of Multidimensional Functional Data objects
#'
#' @description
#' Constructor for `mfd` objects (same as Mfd(...) )
#' 
#' @param argval A list of numeric vectors of argument values at which the `mfd` object is to be evaluated
#' @param X A numeric matrix corresponds to basis expansion coefficients
#'        if `method="coefs"` and discrete observations if `method="data"`.
#' @param mdbs a basismfd object
#' @param method determine the `X` matrix type as "coefs" and "data".
#' @export
Mfd <- function(argval = NULL, X, mdbs, method = "data") mfd$new(argval, X, mdbs, method)
#' @rdname mfd
#' @seealso \code{\link{basismfd}}
