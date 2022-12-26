fpca <- function(fd,method,alpha_orth,ncomp,lambda = NULL,center,lambda_type = "variable"){
  result = ReMFPCA:::fpca_class$new(fd,method,alpha_orth,ncomp,lambda,center,lambda_type = lambda_type)
  #return(fpca_class$new(fd,method,alpha_orth,ncomp,lambda,center))
  return(result)
}