#' @export
"[.mvmfd" <- function(mvmfd_obj, i = "index") {
  if(max(i)> mvmfd_obj$nvar | min(i)<1) stop(" subscript i out of bounds")
  # bs <- mvfd1$basis$basis[[1]].clone()
  
  if(length(i)==1){
    bs <- mvmfd_obj$basis$basis[[i]]
    coef <- mvmfd_obj$coefs[[i]]
    return(mfd$new(X=coef,mdbs=bs,method="coefs"))
  }else{
    bs <- mvmfd_obj$basis$basis[i]
    coef <- mvmfd_obj$coefs[i]
    mvmfd_list <- list()
    for(j in 1:length(i)){
      mvmfd_list[[j]] <- mfd$new(X=coef[[j]],mdbs=bs[[j]],method="coefs")
    }
    return(mvmfd$new(mvmfd_list))
  }
}



#' @export
"[.mvbasismfd" <- function(mvbasismfd_obj, i = "index") {
  if(max(i)> mvbasismfd_obj$nvar | min(i)<1) stop(" subscript i out of bounds")
  if(length(i)==1){
    return(mvbasismfd_obj$basis[[i]])
  }else{
    mvbasismfd_list <- list()
    for(j in 1:length(i)) {
      mvbasismfd_list[[j]] <- mvbasismfd_obj$basis[[j]]
    }    
    return(mvbasismfd$new(mvbasismfd_list))
  }
}
