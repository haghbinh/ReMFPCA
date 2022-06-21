#' A Class of Multidimensional Functional Data objects
#' @description The \code{basismfd} class represents functional data ...
#' @field basis a list of basisfd objects
#' @field dimSupp  a positive integer specify the dimension...
#' 
#' @examples
#' x <- 1
#' 
#' @importFrom fda is.basis eval.basis
#' 
#' @export
mvbasismfd <- R6::R6Class("mvbasismfd",
                      public = list(
                        basis = NULL,
                        dimSupp = NULL,
                        nbasis = list(),
                        supp = list(),
                        nvar = NULL,
                        #' @description
                        #' Constructor for mvbasismfd objects
                        #' @param basis a list of basisfd objects
                        initialize = function(basis) {
                          init_mvbasismfd_check(basis)
                          if (is.basis(basis)) basis <- list(basis)
                          p <- length(basis)
                          for (i in 1:p) {
                            if(is.basis(basis[[i]])){
                              basis[[i]] <- basismfd$new(basis[[i]])
                            }
                            self$dimSupp[i] <- basis[[i]]$dimSupp
                            self$nbasis[[i]] <- basis[[i]]$nbasis
                            self$supp[[i]] <- basis[[i]]$supp 
                          }
                          self$basis <- basis
                          self$nvar <- p
                        },
                        #' @description evalbasismfd
                        #' @param evalarg a list of numeric vector of argument values at which the \code{basismfd} is to be evaluated.
                        #' @return a list
                        eval = function(evalarg) {
                          eval_mvbasismf_validity_check(evalarg,self$nvar)
                          if(is.numeric(evalarg)){
                            evalarg <- list(list(evalarg))
                          }
                          out <- list()
                          for (i in 1:length(evalarg)) {
                            out[[i]]  <- (self$basis[[i]])$eval(evalarg[[i]])
                          }
                          return(out)
                        }
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
    if(self$nvar != 1){
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