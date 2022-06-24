#' @export
length.mfd <-  function(mfd_obj) {
  return(mfd_obj$nobs)
}

#' @export
plot.mfd <- function(mfd_obj,obs=1,xlab="",ylab="",main="",...){
  dimSupp <- mfd_obj$basis$dimSupp
  supp <- mfd_obj$basis$supp
  x_grids <- seq(supp[1,1],supp[2,1],len = 100)
  if(dimSupp == 1){
    X <- mfd_obj$eval(x_grids)
    matplot(x_grids,X,type="l",  lty=1,xlab=xlab,ylab=ylab,...)
  }else{
    y_grids <- seq(supp[1,2],supp[2,2],len = 100)
    X <- mfd_obj$eval(list(x_grids,y_grids))[,,obs]
    image(X ,xlab=xlab,ylab=ylab , axes = FALSE, main=paste(main," Observation:",obs))
    axis(side = 1,at = seq(0,1,len=10),  labels = round(seq(supp[1,1],supp[2,1],len = 10),1))
    axis(side = 2,at = seq(0,1,len=10),  labels =  round(seq(supp[1,2],supp[2,2],len = 10),1))
  }
}



#' @export
"+.mfd" <- function(obj1,obj2) {
  if(length(obj1) != length(obj2)) stop("Two objects must have same length")
  if( obj1$basis$dimSupp != obj2$basis$dimSupp ) stop("Two objects must have same basis dimSupp.")
  if(obj1$basis$dimSupp == 1) {
    if( obj1$basis$basis[[1]]$type != obj2$basis$basis[[1]]$type ) stop("Two objects must have same basis types.")
    bs <- basismfd$new(obj1$basis$basis[[1]])
  }else{
    bs <- list()
    for (i in 1:obj1$basis$dimSupp) {
      if( obj1$basis$basis[[i]]$type != obj2$basis$basis[[i]]$type ) stop("Two objects must have same basis types.")
      bs[[i]] <- basismfd$new(obj1$basis$basis[[i]])
    }
  }
  coef <- obj1$coefs+obj2$coefs  
  return(mfd$new(X=coef,mdbs=bs,method="coefs")) 
}


