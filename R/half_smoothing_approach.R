half_smoothing_approach <- function(data,lambda = NULL,n,center){
  library(expm)
  m.rep = data$nobs
  p = data$nvar
  pc = list()
  #lsv = c()
  #sigma = vector()
  B =c()
  if (is.null(lambda)) {
    for (i in 1:p) {
      lambda = c(lambda,list(2^seq(-20,20,length.out=5)))
    }
  }
  if (p==2) {
    gcv_row = length(lambda[[1]])
    gcv_column = length(lambda[[2]])
  }
  lambda = expand.grid(lambda)
  if (center == TRUE) {
    for (i in 1:p) {
      if (is.matrix(data$coefs[[i]])) {
        c = data$coefs[[i]]-rowMeans(data$coefs[[i]])
        B = rbind(B,c)
      }
      else{
        cc = apply(data$coefs[[i]],3,as.vector)
        c = cc-rowMeans(cc)
        B = rbind(B,c)
      }
    }
  }
  else{
    for (i in 1:p) {
      if (is.matrix(data$coefs[[i]])) {
        c = data$coefs[[i]]
        B = rbind(B,c)
      }
      else{
        cc = apply(data$coefs[[i]],3,as.vector)
        B = rbind(B,cc)
      }
    }
  }
  B = t(as.matrix(B))
  penalty = ReMFPCA:::pen_fun(data,type = "basispen")
  G = data$basis$gram
  G_half = sqrtm(G)
  G_half_inverse = solve(G_half)
  GCV_score = 10^60
  GCV_result = list()
  B_subtilde = B%*%G_half
  
  GCVs = c()
  
  for (j in 1:dim(lambda)[1]) {
    print(j)
    I = ReMFPCA:::I_alpha(data,lambda[j,])
    D = I%*%penalty
    #S = G_half%*%solve(G+D)%*%G_half
    S = sqrtm(solve(G+D))
    #S_half = sqrtm(S)
    #B_tilde = B%*%G_half%*%S_half%*%solve(G_half)
    B_tilde = as.matrix(B%*%G%*%S%*%G_half_inverse)
    X = G_half%*%t(B_tilde)
    #start_time <- Sys.time()
    u=svd(X)$u[,1:n]
    v_temp = svd(X)$v[,1:n]%*%diag(svd(X)$d[1:n])
    b_tilde = G_half_inverse%*%u
    b_temp = S%*%G_half%*%b_tilde
    #GCV = ReMFPCA:::GCVs_rc(data,G,B_tilde_gcv,penalty,p,v_temp,lambda[j,])
    s_alpha_tilde = G_half%*%(S%*%S)%*%G_half
    
    if (all(lambda[j,] == 0)) {
      GCV_score_temp = 0
    }
    else{
      GCV_score_temp = (sum(((diag(dim(s_alpha_tilde)[1])-s_alpha_tilde)%*%(t(B_subtilde)%*%v_temp))^2)/((1-sum(diag(s_alpha_tilde))/dim(G)[1])^2))/dim(G)[1]
    }
    
    GCVs = c(GCVs,GCV_score_temp)
    if (GCV_score_temp < GCV_score) {
      b = b_temp
      v = v_temp
      GCV_score = GCV_score_temp
      GCV_result = lambda[j,]
    }
  }
  #GCVs = log(GCVs)
  if (p == 2) {
    GCVs = matrix(GCVs,nrow = gcv_row,ncol = gcv_column)
  }
  temp_count = 0
  for (j in 1:p) {
    index_start = (temp_count+1)
    index_end = (temp_count+prod(data$basis$nbasis[[j]]))
    pc[[j]] = b[index_start:index_end,]
    temp_count = index_end
  }
  temp = as.matrix(B%*%G%*%b)
  #sigma = sqrt(diag(t(B%*%G%*%b)%*%(B%*%G%*%b)))/sqrt(data$nobs-1)
  sigma = sqrt(diag(t(temp)%*%(temp)))/sqrt(data$nobs-1)
  lsv = (temp)%*%solve(diag(sqrt(diag(t(temp)%*%(temp)))))
  bbbb=c()
  for (kk in 1:p) {
    bbbb = rbind(bbbb,pc[[kk]])
  }
  print(t(bbbb)%*%(G+ReMFPCA:::I_alpha(data,GCV_result)%*%penalty)%*%bbbb)
  return(list(pc,lsv,sigma,GCV_result,GCVs))
}
