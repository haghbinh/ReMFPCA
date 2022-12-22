#' @importFrom Matrix bdiag

I_alpha <- function(alpha,d){
  p <- length(alpha)
  S <- alpha[1]*diag(d[1])
  if (p>1)
    for (i in 2:p) 
      S <- bdiag(S,alpha[i]*diag(d[i]))
  return(S)
}

B_coef <- function(mvmfd_obj){
  p <- mvmfd_obj$nvar
  N <- mvmfd_obj$nobs
  m <- sum(unlist(mvmfd_obj$basis$nbasis))
  B <- matrix(NA, nrow = N, ncol = m)
  for (i in 1:N) {
    mfd_i <- mvmfd_obj[i,]
    B[i,] <- unlist(mfd_i$coefs)
  }
  return(B)
}


slice<-function(X,index) {
  p <- length(index)
  index <- c(0,index)
  n<-nrow(X)
  lapply(seq(1,p),function(i) X[(index[i]+1):min(index[i+1],n),])
}
#conbert Remfpca coefs to mvmfd object
coefs2mvmfd <- function(coefs,mbsfd){
  p <- length(mbsfd$basis)
  n <- ncol(coefs)
  nbasis <- unlist(mbsfd$nbasis)
  X <- slice(coefs,cumsum(nbasis))
  mvlist <- list()
  for (i in 1:p) {
    bs <- mbsfd$basis[[i]]
    mvlist[[i]] <- mfd$new(X=X[[i]], mdbs =bs,method = 'coefs')
  }
  if(p==1){
    return(mvlist[[1]])
  }else{
  return(Mvmfd(mvlist))
  }
}

mvcenteriezed <- function(mvmfd_obj){
  p <- mvmfd_obj$nvar
  cof <- lapply(1:p,function(j) t(scale(t(mvmfd_obj$coefs[[j]]),scale = F)))
  mvlist <- list()
  for (i in 1:p) {
    bs <- mvmfd_obj$basis$basis[[i]]
    mvlist[[i]] <- mfd$new(X=cof[[i]], mdbs =bs,method = 'coefs')
  }
  return(Mvmfd(mvlist))
}


