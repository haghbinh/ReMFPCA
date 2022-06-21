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
basismfd <- R6::R6Class("basismfd",
                      public = list(
                        basis = NULL,
                        dimSupp = NULL,
                        nbasis = NULL,
                        supp = NULL,
                        #' @description
                        #' Constructor for basismfd objects
                        #' @param basis a list of basisfd objects
                        initialize = function(basis) {
                          if(!is.basis(basis) & is.list(basis)) if(!all(sapply(basis, is.basis))){
                            stop("All the elements of basis list must be basisfd object.")
                          }
                          if(is.basis(basis)){
                            self$basis <- list(basis)
                            self$dimSupp <- 1
                            self$supp <- matrix(basis$rangeval,nrow=2,ncol=1)
                            self$nbasis <- basis$nbasis
                          }else{
                            self$basis <- basis
                            self$dimSupp <- length(basis)
                            self$supp <- matrix(0,nrow=2,ncol=length(basis))
                            for (i in 1:length(basis)) {
                              self$nbasis[i] <- basis[[i]]$nbasis
                              self$supp[,i] <- basis[[i]]$rangeval
                            }
                          }
                        },
                        #' @description evalbasismfd
                        #' @param evalarg a list of numeric vector of argument values at which the \code{basismfd} is to be evaluated.
                        #' @return a list
                        eval = function(evalarg) {
                          if(!is.list(evalarg)& !is.numeric(evalarg)){
                            stop('evalarg must be a list or numeric vector')
                          }
                          if(is.numeric(evalarg)){
                            evalarg <- list(evalarg)
                          }
                          if(!all(sapply(evalarg, is.numeric))){
                            stop('evalarg must be a list of numeric vector')
                          }
                          if(length(evalarg) != self$dimSupp) {
                              stop(' length of evalarg list must be equal to dimSupp')
                            }
                          out <- list()
                          for (i in 1:length(evalarg)) {
                            out[[i]] <- eval.basis(evalarg[[i]],self$basis[[i]])
                          }
                          return(out)
                          }
                      )
)