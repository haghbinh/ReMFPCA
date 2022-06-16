# Utility functions

onedG <- function(A, B, grid, method = "trapezoidal") {
  M <- nrow(grid)
  dx <- grid[2] - grid[1]
  if (method == "rectangular") {
    G <- dx * t(B) %*% A
  } else {
    D_diag <- c(1, rep(x = 2, (M - 2)), 1)
    D <- matrix(data = 0, nrow = M, ncol = M)
    diag(D) <- D_diag
    G <- (dx / 2) * t(B) %*% D %*% A
  }
  return(G)
}

twodG <- function(A, B, grid, method = "trapezoidal") {
  if (grid[2, 1] == grid[1, 1]) { # case for by columns
    dy <- grid[2, 2] - grid[1, 2]
    M_1 <- length(which(grid[, 1] == grid[1, 1]))
    dx <- grid[(M_1 + 1), 1] - grid[1, 1]
    M_2 <- length(which(grid[, 2] == grid[1, 2]))
  } else { # case for by rows
    dx <- grid[2, 1] - grid[1, 1]
    M_2 <- length(which(grid[, 2] == grid[1, 2]))
    dy <- grid[(M_2 + 1), 2] - grid[1, 2]
    M_1 <- length(which(grid[, 1] == grid[1, 1]))
  }
  M <- M_1 * M_2
  if (method == "rectangular") {
    G <- (dx * dy) * t(B) %*% A
  } else {
    a_1 <- c(1, rep(x = 2, (M_1 - 2)), 1)
    a_2 <- c(2, rep(x = 4, (M_1 - 2)), 2)
    D_diag <- c(a_1, rep(a_2, (M_2 - 2)), a_1)
    D <- matrix(data = 0, nrow = M, ncol = M)
    diag(D) <- D_diag
    G <- (dx * dy / 4) * t(B) %*% D %*% A
  }
  
  return(G)
}

gram_matrix = function(data,method = "trapezoidal"){
  p = length(data@C)
  G = c()
  for(i in 1:2) {
    #basis functions of functions in different domain
    V = data@B[[i]]
    # 1 dimensional case
    if (dim(data@grid[[i]])[2] == 1) {
      G_sub = ReFPCA:::onedG(A = data@B[[i]],B = data@B[[i]],grid = data@grid[[i]],method = method)
    }
    # 2 dimensional case
    else if (dim(data@grid[[i]])[2] == 2) {
      G_sub = ReFPCA:::twodG(A = data@B[[i]],B = data@B[[i]],grid = data@grid[[i]],method = method)
    }
    #generate block diagonal matrix
    if (i == 1) {
      G = G_sub
    }
    else{
      G = rbind(cbind(G,matrix(0,nrow = nrow(G),ncol = ncol(G_sub))),cbind(matrix(0,nrow = nrow(G_sub),ncol = ncol(G)),G_sub))
    }
  }
  return(G)
}

B_matrix = function(data){
  p = length(data@C)
  B = c()
  for(i in 1:p) {
    #dimean coefficients of basis functions in different domain 
    c = data@C[[i]]-rowMeans(data@C[[i]])
    #generate coefficient matrix 
    B = rbind(B,c)
  }
  B = t(B)
  return(B)
}

I_alpha = function(data,lambda){
  p = length(data@C)
  I = c()
  for(i in 1:p) {
    #I is the diagonal matrix with penalty
    I_sub = lambda[i][[1]]*diag(dim(data@B[[i]])[2])
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

GCVs_rc = function(data,G,B_tilde,penalty,p,v,lambda = NULL){
  gcv = 10^60
  gcv_data = c()
  G_half = eigen(G)$vectors%*%diag(sqrt(eigen(G)$values))%*%t(eigen(G)$vectors)
  if (is.null(lambda)) {
    for (i in 1:p) {
      lambda = c(lambda,list(2^seq(-20,20,length.out=40)))
    }
    lambda = expand.grid(lambda)
    for (j in 1:dim(lambda)[1]) {
      s_alpha = G_half%*%solve(G+ReFPCA:::I_alpha(data,lambda[j,])%*%penalty)%*%G_half
      gcv_temp = (((norm((diag(dim(s_alpha)[1])-s_alpha)%*%(t(B_tilde)%*%v),type = "2"))^2)/((1-sum(diag(s_alpha))/dim(G)[1])^2))/dim(G)[1]
      gcv_data = append(gcv_data,gcv_temp)
      if (gcv_temp < gcv) {
        gcv = gcv_temp
        lambda_target = lambda[j,]
      }
    }
    gcv_data = matrix(gcv_data,40,40)    
    #persp(x = seq(0,40,length.out=nrow(gcv_data)),y = seq(0,40,length.out=ncol(gcv_data)),t(log(gcv_data)))
    return(list(lambda_target,gcv_data))
  }
  else {
    if (length(lambda)!=p) {
      return("wrong input")
    }
    else{
      lambda = expand.grid(lambda)
      for (j in 1:dim(lambda)[1]) {
        s_alpha = G_half%*%solve(G+ReFPCA:::I_alpha(data,lambda[j,])%*%penalty)%*%G_half
        gcv_temp = (((norm((diag(dim(s_alpha)[1])-s_alpha)%*%(t(B_tilde)%*%v),type = "2"))^2)/((1-sum(diag(s_alpha))/dim(G)[1])^2))/dim(G)[1]
        gcv_data = append(gcv_data,gcv_temp)
        if (gcv_temp < gcv) {
          gcv = gcv_temp
          lambda_target = lambda[j,]
        }
      }
      gcv_data = matrix(gcv_data,40,40)
      #persp(x = seq(0,40,length.out=nrow(gcv_data)),y = seq(0,40,length.out=ncol(gcv_data)),t(log(gcv_data)))
      return(list(lambda_target,gcv_data))
    }
  }
}







