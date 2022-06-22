#' A Class of Multivariate Multidimensional Basis Functions
#' @description The \code{mvbasismfd} class represents functional data ...
#' @field basis a list of mvbasisfd objects
#' @field dimSupp  a sequence of positive integer specify the dimension...
#' 
#' @examples
#' x <- 1
#' 
#' @importFrom fda is.basis eval.basis
#' 
#' @export
mvbasismfd <- R6::R6Class("mvbasismfd",
                      public = list(
                        #' @description
                        #' Constructor for mvbasismfd objects
                        #' @param basis a list of basisfd objects
                        initialize = function(basis) {
                          init_mvbasismfd_check(basis)
                          if (is.basis(basis)) basis <- list(basis)
                          private$.nvar <- length(basis)
                          for (i in 1:private$.nvar) {
                            if(is.basis(basis[[i]])){
                              basis[[i]] <- basismfd$new(basis[[i]])
                            } else {
                              basis[[i]] <- basis[[i]]$clone()
                            }
                            private$.dimSupp[i] <- basis[[i]]$dimSupp
                            private$.nbasis[[i]] <- basis[[i]]$nbasis
                            private$.supp[[i]] <- basis[[i]]$supp 
                          }
                          private$.basis <- basis
                        },
                        #' @description evalbasismfd
                        #' @param evalarg a list of numeric vector of argument values at which the \code{basismfd} is to be evaluated.
                        #' @return a list
                        eval = function(evalarg) {
                          eval_mvbasismf_validity_check(evalarg,private$.nvar)
                          if(is.numeric(evalarg)){
                            evalarg <- list(list(evalarg))
                          }
                          out <- list()
                          for (i in 1:length(evalarg)) {
                            out[[i]]  <- (private$.basis[[i]])$eval(evalarg[[i]])
                          }
                          return(out)
                        }
                      ),
                      active = list(
                        nvar = function(value) {
                          if (missing(value)) {
                            private$.nvar
                          } else {
                            stop("`$nvar` is read only", call. = FALSE)
                          }
                        },
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
                        }
                      ),
                      private = list(
                        .nvar = NULL,
                        .basis = NULL,
                        .dimSupp = NULL,
                        .nbasis = list(),
                        .supp = list()
                      )
)



# a function to check the validity of initializer
init_mvbasismfd_check <- function(basis){
  if(is.basis(basis)) {
    basis <- list(basismfd$new(basis))
  }
  if(inherits(basis,"basismfd")) {
    basis <- list(basis)
  }  
  if(is.list(basis)) if(!all(sapply(basis, function(x){
    return(is.basis(x)|inherits(x,"basismfd"))
  }))){
    stop("All the elements of basis list must be basisfd or basismfd object.")
  }
}



# a function to check the validity of evaluation
eval_mvbasismf_validity_check <- function(evalarg,nvar){
  if(!is.list(evalarg)& !is.numeric(evalarg)){
    stop('evalarg must be a list or numeric vector')
  }
  if(is.numeric(evalarg)){
    if(nvar != 1){
      stop('evalarg is allowd be a numeric if nvar = 1.')
    }else{
      evalarg <- list(list(evalarg))
    }
  }
  if(!all(sapply(evalarg, function(x) is.numeric(x) | is.list(x)))){
    stop('evalarg list elements must be a list or numeric vector')
  }
  if(length(evalarg) != nvar){
    stop('length of evalarg is not equal to nvar.')
  }
}