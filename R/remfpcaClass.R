#' A Class for ReMFPCA objects
#' @description The \code{remfpca} class represents functional data ...
#' @field basis a basismfd object
#' @field coeff a matrix with nrow=subjects and ncol=total number of basis ...
#' 
#' @examples
#' x <- 1
#' 
#' @importFrom fda getbasispenalty
#' @importFrom Matrix bdiag
#' @importFrom expm sqrtm 
#' @export
remfpca <- R6::R6Class("remfpca",
    public = list(
      #' @description
      #' Constructor for mfd objects
      #' @param mdbs a basismfd object
      initialize = function(mvmfd_obj, method="eigen", alpha=NULL, ncomp=2, centerfns = FALSE) {
           if (is.mfd(mvmfd_obj)) mvmfd_obj <- mvmfd$new(mvmfd_obj)
           p <- mvmfd_obj$nvar
           N <- mvmfd_obj$nobs
           nbasis <- unlist(mvmfd_obj$basis$nbasis)
           if (is.null(alpha)) 
             alpha <- rep(0,len=p)
           if(centerfns == TRUE){
             mvmfd_obj <- mvcenteriezed(mvmfd_obj)
           }
           d <- c(NA)
           G <- D <- list()
           for (i in 1:p) {
             mfd_i <- mvmfd_obj[,i]
             d[i] <- nrow(mfd_i$coefs)
             G[[i]] <- mfd_i$basis$gram
             D[[i]] <- getbasispenalty(mfd_i$basis$basis[[1]])
             }
           I_a <- I_alpha(alpha,d)
           G <- as.matrix(bdiag(G))
           D <- bdiag(D)
           # E <- eigen(G+I_a%*%D)
           # S <- solve((E$vectors)%*%diag(sqrt(E$values))%*%t(E$vectors))
           S <- solve(sqrtm(G+I_a%*%D))
           B <- B_coef(mvmfd_obj)
           V <- 1/(N-1)*t(B)%*%(diag(1,N)-1/N*matrix(1,nrow = N, ncol = N))%*%B
           Eig <- eigen(S%*%t(G)%*%V%*%G%*%t(S))
           coefs <- t(S)%*%(Eig$vectors)[,1:ncomp]
           private$.values <- Eig$values[1:ncomp]
           mbsfd <- mvmfd_obj$basis
           private$.Eigenfunctions <- coefs2mvmfd(coefs,mbsfd)
           }
      ),
    active = list(
      Eigenfunctions = function(value) {
        if (missing(value)) {
          private$.Eigenfunctions
        } else {
          stop("`$Eigenfunctions` is read only", call. = FALSE)
        }
      },
      values = function(value) {
        if (missing(value)) {
          private$.values
        } else {
          stop("`$values` is read only", call. = FALSE)
        }
      }
    ),
    
    private = list(
      .Eigenfunctions = NULL,
      .values = NULL
      )
)

#' @export
Remfpca <- function(mvmfd_obj, method="eigen", alpha=NULL, ncomp=2, center = TRUE) remfpca$new(mvmfd_obj, method, alpha, ncomp, center)
