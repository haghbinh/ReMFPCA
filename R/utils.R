I_alpha <- function(data, lambda) {
  p <- length(data$coefs)
  I <- c()
  for (i in 1:p) {
    # I is the diagonal matrix with penalty
    I_sub <- lambda[i][[1]] * diag(prod(data$basis$nbasis[[i]]))
    # generate block diagonal matrix
    if (i == 1) {
      I <- I_sub
    } else {
      I <- rbind(cbind(I, matrix(0, nrow = nrow(I), ncol = ncol(I_sub))), cbind(matrix(0, nrow = nrow(I_sub), ncol = ncol(I)), I_sub))
    }
  }
  return(I)
}


# GCVs_rc <- function(data, G, G_half, B_tilde, penalty, p, v, lambda) {
#   s_alpha_tilde <- G_half %*% solve(G + I_alpha(data, lambda) %*% penalty) %*% G_half
#   gcv_score <- (sum(((diag(dim(s_alpha_tilde)[1]) - s_alpha_tilde) %*% (t(B_tilde) %*% v))^2) / ((1 - sum(diag(s_alpha_tilde)) / dim(G)[1])^2)) / dim(G)[1]
#   return(gcv_score)
# }
