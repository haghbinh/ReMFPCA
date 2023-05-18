#' @title A Class of Multidimensional Functional Data objects
#'
#' @description
#' The 'mvmfd' class represents functional data ...
#'
#' @field basis A `basismfd` object
#' @field coefs A matrix with nrow = subjects and ncol = total number of basis ...
#' @field nobs description
#' @field nvar description
#' 
#' @examples
#' require(fda)
#' bs1 <- create.fourier.basis(c(0, 2 * pi), 5)
#' bs2 <- create.bspline.basis(c(0, 1), 7)
#' bs3 <- create.exponential.basis(c(0, 2), 3)
#' argval1 <- seq(0,2*pi,len=12)
#' X1 <- outer(sin(argval1),seq(0.5,1.5,len=nobs))
#' mdbs1 <- Basismfd(bs1)
#' mfd1 <- Mfd(argval1,X1,mdbs1)
#' mdbs2 <- Basismfd(bs1)
#' argval2 <- argval1
#' X2 <- outer(cos(argval2),seq(0.2,1.5,len=nobs))
#' mfd2 <- Mfd(argval2,X2,mdbs1)
#' mvmfd1 <- Mvmfd(mfd1,mfd2)
#' mvmfd1[1]
#' mvmfd1[1,1]
#' mvmfd1[1:5,2]
#' mvmfd1[,1]
#' mvmfd1[1:5,]
#' evalarg <- list(argval1,argval2)
#' mvmfd1$eval(evalarg)
#' mvmfd1+mvmfd1
#' mean(mvmfd1)
#' inprod_mvmfd(mvmfd1,mvmfd1)
#' norm_mvmfd(mvmfd1)
#' plot(mvmfd1)
#' bimfdplot(mvmfd1)
#' 
#' @seealso \code{\link{basismfd}} 
#' 
#' @importFrom fda is.basis eval.basis Data2fd
#'
#' @export
mvmfd <- R6::R6Class("mvmfd",
                     public = list(
                       #' @param ... A list of `mfd` objects
                       #'
                       initialize = function(...) {
                         mfd_list <- list(...)
                         if (is.list(mfd_list[[1]])) mfd_list <- mfd_list[[1]]
                         # if (is.mfd(mfd_list)) mfd_list <- list(mfd_list)
                         init_mvmfd_check(mfd_list)
                         basis_list <- list()
                         private$.nobs <- mfd_list[[1]]$nobs
                         private$.nvar <- length(mfd_list)
                         for (i in 1:private$.nvar) {
                           mfd_list[[i]] <- mfd_list[[i]]$clone()
                           basis_list[[i]] <- mfd_list[[i]]$basis
                           private$.coefs[[i]] <- mfd_list[[i]]$coefs
                         }
                         private$.basis <- mvbasismfd$new(basis_list)
                       },
                       
                       #' @description
                       #' Eval method for `mvmfd` objects
                       #'
                       #' @param evalarg A list of numeric vectors of argument values at which the `mvmfd` is to be evaluated.
                       #' @return A list of evaluated values
                       eval = function(evalarg) {
                         eval_mvmfd_validity_check(evalarg)
                         Bmat <- private$.basis$eval(evalarg)
                         Xhat <- list()
                         for (i in 1:private$.nvar) {
                           if (is.matrix(private$.coefs[[i]])) {
                             Xhat[[i]] <- Bmat[[i]][[1]] %*% private$.coefs[[i]]
                           } else {
                             Xhat[[i]] <- (Bmat[[i]][[2]] %x% Bmat[[i]][[1]]) %*% apply(private$.coefs[[i]], 3, as.vector)
                             Xhat[[i]] <- array(Xhat[[i]], dim = c(sapply(evalarg[[i]], length), private$.nobs))
                           }
                         }
                         return(Xhat)
                       },
                       
                       #' @description
                       #' Print method for `mvmfd` objects
                       #'
                       #' @param ... Additional arguments to be passed to `print`
                       #'
                       print = function(...) {
                         cat("A 'mvmfd' object with", private$.nvar, "variable(s):\n")
                         for (i in 1:private$.nvar) {
                           cat("\nVariable ", i, ":\n", sep = "")
                           print(self[, i])
                         }
                         invisible(self)
                       }
                     ),
                     active = list(
                       # basis field
                       basis = function(value) {
                         if (missing(value)) {
                           private$.basis
                         } else {
                           stop("`$basis` is read only", call. = FALSE)
                         }
                       },
                       
                       # coefs field
                       coefs = function(value) {
                         if (missing(value)) {
                           private$.coefs
                         } else {
                           stop("`$coefs` is read only", call. = FALSE)
                         }
                       },
                       
                       # nvar field
                       nvar = function(value) {
                         if (missing(value)) {
                           private$.nvar
                         } else {
                           stop("`$nvar` is read only", call. = FALSE)
                         }
                       },
                       
                       # nobs field
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
                       .coefs = list(),  # we record vectorized of the coefs
                       .nobs = NULL,
                       .nvar = NULL
                     )
)

# A function to check the validity of initializer
init_mvmfd_check <- function(mfd_list) {
  if (!all(sapply(mfd_list, is.mfd))) {
    stop("All the elements of the inputs list must have the class of `mfd`")
  }
  n <- mfd_list[[1]]$nobs
  for (y in mfd_list) {
    if (n != y$nobs) stop("The number of observations in all variables should be equal.")
  }
}

# A function to check the validity of evaluation
eval_mvmfd_validity_check <- function(evalarg, dimSupp) {
  x <- 1
}

#' @export
Mvmfd <- function(...) mvmfd$new(...)
