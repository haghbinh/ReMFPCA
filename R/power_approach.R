#' @importFrom expm sqrtm
power_algo_fun <- function(mvmfd_obj, n, alpha, centerfns, alpha_orth, lambda_type,penalty_type) {
  p <- mvmfd_obj$nvar
  # penalty = pen_fun(mvmfd_obj,type = "basispen")
  penalty <- pen_fun(mvmfd_obj, type = penalty_type)
  B <- c()
  if (is.null(alpha)) {
    for (i in 1:p) {
      alpha <- c(alpha, list(2^seq(-20, 20, length.out = 10)))
    }
  }

  if (lambda_type == "variable") {
    if (p == 2) {
      gcv_row <- length(alpha[[1]])
      gcv_column <- length(alpha[[2]])
    }
    alpha <- expand.grid(alpha)
  } else if (lambda_type == "component") {
    if (p == 2) {
      gcv_row <- 1
      gcv_column <- 1
    }
  } else {
    return("wrong put")
  }

  pc <- list()
  if (centerfns == TRUE) {
    for (i in 1:p) {
      if (is.matrix(mvmfd_obj$coefs[[i]])) {
        c <- mvmfd_obj$coefs[[i]] - rowMeans(mvmfd_obj$coefs[[i]])
        B <- rbind(B, c)
      } else {
        cc <- apply(mvmfd_obj$coefs[[i]], 3, as.vector)
        c <- cc - rowMeans(cc)
        B <- rbind(B, c)
      }
    }
  } else {
    for (i in 1:p) {
      if (is.matrix(mvmfd_obj$coefs[[i]])) {
        c <- mvmfd_obj$coefs[[i]]
        B <- rbind(B, c)
      } else {
        cc <- apply(mvmfd_obj$coefs[[i]], 3, as.vector)
        B <- rbind(B, cc)
      }
    }
  }
  B <- t(B)
  lsv <- c()
  sigma <- vector()
  G <- mvmfd_obj$basis$gram
  G_half <- sqrtm(G)
  if (alpha_orth == FALSE) {
    GCV_socre <- list()
    GCVs <- list()
    I <- list()
    D <- list()
    S <- list()
    s_alpha_tilde <- list()
    for (i in 1:n) {
      GCV_score <- 10^60
      # print("********")
      # print(i)
      # print("********")
      if (i != 1) {
        b <- as.matrix(b)
        B <- as.matrix(B - v %*% t(b))
      }
      B_subtilde <- as.matrix(B %*% G_half)
      b_temp <- svd(B_subtilde)$v[, 1]
      # GCVs[[i]] = c()
      if (!is.null(dim(alpha))) {
        count <- dim(alpha)[1]
      } else {
        count <- 1
      }
      for (j in 1:count) {
        error_b <- 1000
        # print(j)
        # if (i == 1) {
        if (lambda_type == "variable") {
          if (i == 1) {
            I[[j]] <- I_alpha(mvmfd_obj, alpha[j, ])
            D[[j]] <- I[[j]] %*% penalty
            S[[j]] <- sqrtm(solve(G + D[[j]]))
            s_alpha_tilde[[j]] <- G_half %*% (S[[j]] %*% S[[j]]) %*% G_half
          }
        } else {
          I[[i]] <- I_alpha(mvmfd_obj, alpha[[i]])
          D[[i]] <- I[[i]] %*% penalty
          S[[i]] <- sqrtm(solve(G + D[[i]]))
          s_alpha_tilde[[i]] <- G_half %*% (S[[i]] %*% S[[i]]) %*% G_half
        }
        # I[[j]] =  I_alpha(mvmfd_obj,alpha[j,])
        # D[[j]] = I[[j]]%*%penalty
        # S[[j]] = sqrtm(solve(G+D[[j]]))
        # s_alpha_tilde[[j]] = G_half%*%(S[[j]]%*%S[[j]])%*%G_half
        # }
        # I = I_alpha(mvmfd_obj,alpha[j,])
        # D = I%*%penalty
        # S = sqrtm(solve(G+D))
        # s_alpha_tilde = G_half%*%(S%*%S)%*%G_half
        while (error_b > 0.001) {
          if (lambda_type == "variable") {
            b_previous <- b_temp
            v_temp <- (B %*% G %*% b_temp)
            b_temp <- (S[[j]] %*% S[[j]] %*% G %*% t(B) %*% v_temp)
            b_temp <- as.matrix(b_temp)
            b_temp <- b_temp %*% solve(sqrt(t(b_temp) %*% (G + D[[j]]) %*% b_temp))
            error_b <- sum((b_temp - b_previous)^2)
          } else {
            b_previous <- b_temp
            v_temp <- (B %*% G %*% b_temp)
            b_temp <- (S[[i]] %*% S[[j = i]] %*% G %*% t(B) %*% v_temp)
            b_temp <- as.matrix(b_temp)
            b_temp <- b_temp %*% solve(sqrt(t(b_temp) %*% (G + D[[i]]) %*% b_temp))
            error_b <- sum((b_temp - b_previous)^2)
          }
        }

        if (lambda_type == "variable") {
          if (all(alpha[j, ] == 0)) {
            GCV_score_temp <- 0
          } else {
            GCV_score_temp <- (sum(((diag(dim(s_alpha_tilde[[j]])[1]) - s_alpha_tilde[[j]]) %*% (t(B_subtilde) %*% v_temp))^2) / ((1 - sum(diag(s_alpha_tilde[[j]])) / dim(G)[1])^2)) / dim(G)[1]
          }
        } else {
          GCV_score_temp <- (sum(((diag(dim(s_alpha_tilde[[j]])[1]) - s_alpha_tilde[[j]]) %*% (t(B_subtilde) %*% v_temp))^2) / ((1 - sum(diag(s_alpha_tilde[[j]])) / dim(G)[1])^2)) / dim(G)[1]
        }
        if (j == 1) {
          GCVs[[i]] <- GCV_score_temp
        } else {
          GCVs[[i]] <- c(GCVs[[i]], GCV_score_temp)
        }
        if (GCV_score_temp < GCV_score) {
          GCV_score <- GCV_score_temp
          if (lambda_type == "variable") {
            GCV_result <- alpha[j, ]
          } else {
            GCV_result <- alpha[[i]]
          }
          # GCV_result = alpha[j,]
          v <- v_temp
          b <- b_temp
        }
      }
      # GCVs[[i]] = log(GCVs[[i]])
      if (p == 2) {
        GCVs[[i]] <- matrix(GCVs[[i]], nrow = gcv_row, ncol = gcv_column)
      }
      GCV_socre[[i]] <- GCV_result
      temp_count <- 0
      for (j in 1:p) {
        index_start <- (temp_count + 1)
        index_end <- (temp_count + prod(mvmfd_obj$basis$nbasis[[j]]))
        if (i == 1) {
          pc[[j]] <- b[index_start:index_end, ]
        } else {
          pc[[j]] <- cbind(pc[[j]], b[index_start:index_end, ])
        }
        temp_count <- index_end
      }
      lsv_before <- (B %*% G %*% b) / norm(B %*% G %*% b, type = "2")
      lsv <- cbind(lsv, lsv_before)
      sigma[i] <- norm(B %*% G %*% b, type = "2") / sqrt(mvmfd_obj$nobs - 1)
    }

    bbbb <- c()
    for (k in 1:p) {
      bbbb <- rbind(bbbb, pc[[k]])
    }
  } else if (alpha_orth == TRUE) {
    B_subtilde <- as.matrix(B %*% G_half)
    b_temp <- svd(B_subtilde)$v[, 1:n]
    GCVs <- c()
    GCV_score <- 10^60

    for (j in 1:dim(alpha)[1]) {
      # print(j)
      error_b <- 1000
      I <- I_alpha(mvmfd_obj, alpha[j, ])
      D <- I %*% penalty
      S <- sqrtm(solve(G + D))
      s_alpha_tilde <- G_half %*% (S %*% S) %*% G_half
      while (error_b > 0.00001) {
        b_previous <- b_temp
        v_temp <- (B %*% G %*% b_temp)
        b_temp <- S %*% S %*% G_half %*% qr.Q(qr(as.matrix(G_half %*% t(B) %*% v_temp)))
        b_temp <- as.matrix(b_temp)
        temp <- as.matrix((t(b_temp) %*% (G + D) %*% b_temp))
        b_temp <- as.matrix((b_temp) %*% solve(sqrtm(diag(diag(temp)))))
        error_b <- sum((b_temp - b_previous)^2)
      }
      if (all(alpha[j, ] == 0)) {
        GCV_score_temp <- 0
      } else {
        GCV_score_temp <- (sum(((diag(dim(s_alpha_tilde)[1]) - s_alpha_tilde) %*% (t(B_subtilde) %*% v_temp))^2) / ((1 - sum(diag(s_alpha_tilde)) / dim(G)[1])^2)) / dim(G)[1]
      }

      # GCV_score_temp = (sum(((diag(dim(s_alpha_tilde)[1])-s_alpha_tilde)%*%(t(B_subtilde)%*%v_temp))^2)/((1-sum(diag(s_alpha_tilde))/dim(G)[1])^2))/dim(G)[1]
      GCVs <- c(GCVs, GCV_score_temp)
      if (GCV_score_temp < GCV_score) {
        GCV_score <- GCV_score_temp
        GCV_result <- alpha[j, ]
        v <- v_temp
        b <- b_temp
        # print(GCV_score)
      }
    }

    # GCVs = log(GCVs)
    if (p == 2) {
      GCVs <- matrix(GCVs, nrow = gcv_row, ncol = gcv_column)
    }
    GCV_socre <- GCV_result
    temp_count <- 0
    for (j in 1:p) {
      index_start <- (temp_count + 1)
      index_end <- (temp_count + prod(mvmfd_obj$basis$nbasis[[j]]))
      pc[[j]] <- b[index_start:index_end, ]
      temp_count <- index_end
    }
    temp2 <- as.matrix(B %*% G %*% b)
    sigma <- sqrt(diag(t(temp2) %*% (temp2))) / sqrt(mvmfd_obj$nobs - 1)
    lsv <- (temp2) %*% solve(diag(sqrt(diag(t(temp2) %*% (temp2)))))


    bbbb <- c()
    for (k in 1:p) {
      bbbb <- rbind(bbbb, pc[[k]])
    }
    print(t(bbbb) %*% (G + I_alpha(mvmfd_obj, GCV_result) %*% penalty) %*% bbbb)
  }
  return(list(pc, lsv, sigma, GCV_socre, GCVs))
}
