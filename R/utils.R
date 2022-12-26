I_alpha = function(data,lambda){
  p = length(data$coefs)
  I = c()
  for(i in 1:p) {
    #I is the diagonal matrix with penalty
    I_sub = lambda[i][[1]]*diag(prod(data$basis$nbasis[[i]]))
    #generate block diagonal matrix
    if (i == 1) {
      I = I_sub
    }
    else{
      I = rbind(cbind(I,matrix(0,nrow = nrow(I),ncol = ncol(I_sub))),cbind(matrix(0,nrow = nrow(I_sub),ncol = ncol(I)),I_sub))
    }
  }
  return(I)
}

# GCVs_rc = function(data,G,B_tilde,penalty,p,v,lambda = NULL){
#   B_tilde = as.matrix(B_tilde)
#   G = as.matrix(G)
#   gcv = 10^60
#   gcv_data = c()
#   G_half = eigen(G)$vectors%*%diag(sqrt(eigen(G)$values))%*%t(eigen(G)$vectors)
#   if (is.null(lambda)) {
#     for (i in 1:p) {
#       lambda = c(lambda,list(2^seq(-20,20,length.out=40)))
#     }
#     lambda = expand.grid(lambda)
#     for (j in 1:dim(lambda)[1]) {
#       s_alpha_tilde = G_half%*%solve(G+ReMFPCA:::I_alpha(data,lambda[j,])%*%penalty)%*%G_half
#       #gcv_temp = (((norm((diag(dim(s_alpha_tilde)[1])-s_alpha_tilde)%*%(t(B_tilde)%*%v),type = "2"))^2)/((1-sum(diag(s_alpha_tilde))/dim(G)[1])^2))/dim(G)[1]
#       gcv_temp = (sum(((diag(dim(s_alpha_tilde)[1])-s_alpha_tilde)%*%(t(B_tilde)%*%v))^2)/((1-sum(diag(s_alpha_tilde))/dim(G)[1])^2))/dim(G)[1]
#       gcv_data = append(gcv_data,gcv_temp)
#       if (gcv_temp < gcv) {
#         gcv = gcv_temp
#         lambda_target = lambda[j,]
#       }
#     }
#     gcv_data = matrix(gcv_data,40,40)
#     #persp(x = seq(0,40,length.out=nrow(gcv_data)),y = seq(0,40,length.out=ncol(gcv_data)),t(log(gcv_data)))
#     return(list(lambda_target,gcv_data,gcv))
#   }
#   else {
#     if (length(lambda)!=p) {
#       return("wrong input")
#     }
#     else{
#       lambda = expand.grid(lambda)
#       for (j in 1:dim(lambda)[1]) {
#         s_alpha_tilde = G_half%*%solve(G+ReMFPCA:::I_alpha(data,lambda[j,])%*%penalty)%*%G_half
#         #gcv_temp = (((norm((diag(dim(s_alpha_tilde)[1])-s_alpha_tilde)%*%(t(B_tilde)%*%v),type = "2"))^2)/((1-sum(diag(s_alpha_tilde))/dim(G)[1])^2))/dim(G)[1]
#         gcv_temp = (sum(((diag(dim(s_alpha_tilde)[1])-s_alpha_tilde)%*%(t(B_tilde)%*%v))^2)/((1-sum(diag(s_alpha_tilde))/dim(G)[1])^2))/dim(G)[1]
#         gcv_data = append(gcv_data,gcv_temp)
#         if (gcv_temp < gcv) {
#           gcv = gcv_temp
#           lambda_target = lambda[j,]
#         }
#       }
#       gcv_data = matrix(gcv_data,40,40)
#       #persp(x = seq(0,40,length.out=nrow(gcv_data)),y = seq(0,40,length.out=ncol(gcv_data)),t(log(gcv_data)))
#       return(list(lambda_target,gcv_data,gcv))
#     }
#   }
# }

GCVs_rc = function(data,G,B_tilde,penalty,p,v,lambda){
  s_alpha_tilde = G_half%*%solve(G+ReMFPCA:::I_alpha(data,lambda)%*%penalty)%*%G_half
  gcv_score = (sum(((diag(dim(s_alpha_tilde)[1])-s_alpha_tilde)%*%(t(B_tilde)%*%v))^2)/((1-sum(diag(s_alpha_tilde))/dim(G)[1])^2))/dim(G)[1]
  return(gcv_score)
}

