#' @export
length.mfd <-  function(mfd_obj) {
  return(mfd_obj$nobs)
}

#' @export
plot.mfd <- function(mfd_obj,obs=1,xlab="",ylab="",main="",...){
  dimSupp <- mfd_obj$basis$dimSupp
  supp <- mfd_obj$basis$supp
  x_grids <- seq(supp[1,1],supp[2,1],len = 1000)
  if(dimSupp == 1){
    X <- mfd_obj$eval(x_grids)
    matplot(x_grids,X,type="l",  lty=1,xlab=xlab,ylab=ylab,main=main,...)
  }else{
    y_grids <- seq(supp[1,2],supp[2,2],len = 100)
    X <- mfd_obj$eval(list(x_grids,y_grids))[,,obs]
    image(X ,xlab=xlab,ylab=ylab , axes = FALSE, main=paste(main," Observation:",obs))
    axis(side = 1,at = seq(0,1,len=10),  labels = round(seq(supp[1,1],supp[2,1],len = 10),1))
    axis(side = 2,at = seq(0,1,len=10),  labels =  round(seq(supp[1,2],supp[2,2],len = 10),1))
  }
}

#' @export
"+.mfd" <- function(obj1, obj2=NULL) {
  if (is.null(obj2)) return(obj1)
  if (xor(is.mfd(obj1), is.mfd(obj2))) {
    if (!xor(is.double(obj1), is.double(obj2))) stop("At least one object must be an mfd, and the other one can be a scalar")
    if (is.double(obj1)) {
      temp <- obj1
      obj1 <- obj2
      obj2 <- temp
    }
    coef <- obj1$coefs + obj2
  } else {
    if (length(obj1) != length(obj2)) stop("Two objects must have same length")
    if (obj1$basis$dimSupp != obj2$basis$dimSupp) stop("Two objects must have same basis dimSupp.")
    if (obj1$basis$dimSupp == 1) {
      if (obj1$basis$basis[[1]]$type != obj2$basis$basis[[1]]$type) stop("Two objects must have same basis types.")
    } else {
      for (i in 1:obj1$basis$dimSupp) {
        if (obj1$basis$basis[[i]]$type != obj2$basis$basis[[i]]$type) stop("Two objects must have same basis types.")
      }
    }
    coef <- obj1$coefs + obj2$coefs
  }
  return(mfd$new(X = coef, mdbs = obj1$basis$clone(), method = "coefs"))
}

#' @export
"*.mfd" <- function(obj1, obj2) {
  if (xor(is.mfd(obj1), is.mfd(obj2))) {
    if (xor(is.double(obj1), is.double(obj2)))
    if (is.double(obj1)) {
      temp <- obj1
      obj1 <- obj2
      obj2 <- temp
    }
    coef <- obj1$coefs * obj2
  } else 
    stop("One object must be an mfd, and the other one a scalar")
  return(mfd$new(X = coef, mdbs = obj1$basis$clone(), method = "coefs"))
}


#' @export
"-.mfd" <- function(obj1, obj2=NULL) {
  if (is.null(obj2)) return((-1)*obj1)
  if (xor(is.mfd(obj1), is.mfd(obj2))) {
    if (!xor(is.double(obj1), is.double(obj2))) stop("At least one object must be an mfd, and the other one can be a scalar")
    if (is.double(obj1)) {
      coef <- obj1 + (-1)*obj2$coefs
      obj1 <- obj2
    } else{
      coef <- obj1$coefs + (-obj2)
    }
  } else return(obj1 + (-1)*obj2)
  return(mfd$new(X = coef, mdbs = obj1$basis$clone(), method = "coefs"))
}

#' @export
mean.mfd <- function(mfd_obj){
  cof <- apply(mfd_obj$coefs,2,mean)
  bs <- mfd_obj$basis
  return(Mfd(X=cof,mdbs=bs,method = 'coefs'))
}


#' @export
mean.mvmfd <- function(mvmfd_obj){
  p <- mvmfd_obj$nvar
  mvlist <- lapply(1:p,function(j) mean(mvmfd_obj[,j]))
  return(Mvmfd(mvlist))
}

# inprod <- function(object, ...) UseMethod("inprod")
#' @importFrom fda fd inprod
#' @export
inprod_mfd <- function(mfd_obj1,mfd_obj2){
  bs1 <- mfd_obj1$basis$basis[[1]]
  bs2 <- mfd_obj2$basis$basis[[1]]
  cof1 <- mfd_obj1$coefs
  cof2 <- mfd_obj2$coefs
  fd1 <- fd(coef = cof1,basisobj = bs1)
  fd2 <- fd(coef = cof2,basisobj = bs2)
  inpr <- fda::inprod(fd1,fd2)
  return(inpr)
}
  

#' @export
inprod_mvmfd <- function(mvmfd_obj1,mvmfd_obj2){
  p <- mvmfd_obj1$nvar
  if (p!= mvmfd_obj2$nvar) stop('The number of variablles must be equal.')
  m <- mvmfd_obj1$nobs
  n <- mvmfd_obj2$nobs
  inpr <- matrix(0,nrow = m,ncol = n)
  for (j in 1:p) {
    inpr <- inpr+inprod_mfd(mvmfd_obj1[,j],mvmfd_obj2[,j])
  }
  return(inpr)
}


