#' @export
"[.mfd" <- function(mfd_obj, i = "index") {
  if(is.null(i)) i <- 1:mfd_obj$nobs
  if (max(i) > mfd_obj$nobs | min(i) < 1) stop(" subscript i out of bounds")
  bs <- mfd_obj$basis$clone()
  if(mfd_obj$basis$dimSupp==1){
    coef <- mfd_obj$coefs[, i]
  }else{
    coef <- mfd_obj$coefs[, , i]
  }
  return(mfd$new(X = coef, mdbs = bs, method = "coefs"))
}



#' @export
"[.mvmfd" <- function(mvmfd_obj, i = "index", j = "index") {
  if (i[1] == "index") i <- 1:mvmfd_obj$nobs
  if (j[1] == "index") j <- 1:mvmfd_obj$nvar
  if (max(i) > mvmfd_obj$nobs | min(i) < 1) stop(" subscript i out of bounds")
  if (max(j) > mvmfd_obj$nvar | min(i) < 1) stop(" subscript j out of bounds")
  dimSupp <- mvmfd_obj$basis$dimSupp
  bs <- mvmfd_obj$basis[j]
  if (length(j) == 1) {
    if (dimSupp[j] == 1) {
      coef <- mvmfd_obj$coefs[[j]][, i]
    } else {
      coef <- mvmfd_obj$coefs[[j]][, , i]
    }
    return(mfd$new(X = coef, mdbs = bs, method = "coefs"))
  } else {
    mvmfd_list <- list()
    for (k in 1:length(j)) {
      if (dimSupp[k] == 1) {
        coef <- mvmfd_obj$coefs[[k]][, i]
      } else {
        coef <- mvmfd_obj$coefs[[k]][, , i]
      }
      mvmfd_list[[k]] <- mfd$new(X = coef, mdbs = bs[k], method = "coefs")
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
