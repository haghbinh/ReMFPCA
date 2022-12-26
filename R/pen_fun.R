pen_fun = function(data,devorder=2,type){
  D_final = c()
  if (type == 'basispen') {
    for (i in 1:data$nvar) {
      # one dimensional case
      if (data$basis$basis[[i]]$dimSupp == 1) {
        D = fda::eval.penalty(data$basis$basis[[i]]$basis[[1]],Lfdobj = devorder)
      }
      # two dimensional case with fourier or bspline (two basis input)
      else if (data$basis$basis[[i]]$dimSupp == 2) {
        P1 = fda::eval.penalty(data$basis$basis[[i]]$basis[[1]],Lfdobj = devorder)
        P2 = fda::eval.penalty(data$basis$basis[[i]]$basis[[2]],Lfdobj = devorder)
        D = P1%x%diag(nrow(P2))+diag(nrow(P1))%x%P2
      }
      #generate block diagonal matrix
      if (i == 1) {
        D_final = D
      }
      else{
        D_final = rbind(cbind(D_final,matrix(0,nrow = nrow(D_final),ncol = ncol(D))),cbind(matrix(0,nrow = nrow(D),ncol = ncol(D_final)),D))
      }
    }
  }
  else if (type =='coefpen') {
    for (i in 1:data$nvar) {
      # 1-dimensional case
      if (data$basis$basis[[i]]$dimSupp == 1) {
        L = diff(diag(data$basis$basis[[i]]$nbasis),differences=devorder)
        D = t(L)%*%L
      }
      # 2-dimensional case
      else if (data$basis$basis[[i]]$dimSupp == 2) {
        P1 = diff(diag(data$basis$basis[[i]]$nbasis[1]),differences=devorder)
        P2 = diff(diag(data$basis$basis[[i]]$nbasis[2]),differences=devorder)
        D = P1%x%diag(nrow(P2))+diag(nrow(P1))%x%P2
      }
      #generate block diagonal matrix
      if (i == 1) {
        D_final = D
      }
      else{
        D_final = rbind(cbind(D_final,matrix(0,nrow = nrow(D_final),ncol = ncol(D))),cbind(matrix(0,nrow = nrow(D),ncol = ncol(D_final)),D))
      }
    }
  }
  return(D_final)
}
