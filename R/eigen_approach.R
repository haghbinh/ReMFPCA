#' @importFrom expm sqrtm
#' @importFrom utils  txtProgressBar setTxtProgressBar

eigen_approach <-
  function(mvmfd_obj,
           n,
           alpha,
           centerfns,
           penalty_type) {
    m.rep <- mvmfd_obj$nobs
    p <- mvmfd_obj$nvar
    indices <-
      sapply(1:p, function(i)
        prod(mvmfd_obj$basis$nbasis[[i]]))
    if (is.null(alpha)) {
      for (i in 1:p) {
        alpha <- c(alpha, list(2 ^ seq(-20, 20, length.out = 5)))
      }
    }
    if (p == 2) {
      gcv_row <- length(alpha[[1]])
      gcv_column <- length(alpha[[2]])
      index1 = mvmfd_obj$basis$nbasis[[1]]
      index2 = mvmfd_obj$basis$nbasis[[2]]
    }
    alpha <- expand.grid(alpha)
    penalty <- pen_fun(mvmfd_obj, type = penalty_type)
    G <- as.matrix(mvmfd_obj$basis$gram)
    G_half <- expm::sqrtm(G)
    
    B <- c()
    B_c <- c()
    if (centerfns) {
      for (i in 1:p) {
        if (is.matrix(mvmfd_obj$coefs[[i]])) {
          B <- rbind(B, mvmfd_obj$coefs[[i]])
          B_c <-
            rbind(B_c, mvmfd_obj$coefs[[i]] - rowMeans(mvmfd_obj$coefs[[i]]))
        } else {
          cc <- apply(mvmfd_obj$coefs[[i]], 3, as.vector)
          B <- rbind(B, cc)
          B_c <- rbind(B_c, cc - rowMeans(cc))
        }
      }
    } else {
      for (i in 1:p) {
        if (is.matrix(mvmfd_obj$coefs[[i]])) {
          B <- rbind(B, mvmfd_obj$coefs[[i]])
        } else {
          cc <- apply(mvmfd_obj$coefs[[i]], 3, as.vector)
          B <- rbind(B, cc)
        }
      }
      B_c <- B
    }
    B <- t(B)
    B_c <- t(B_c)
    # B_subtilde <- B_c %*% G_half
    I_matrix <- diag(1, m.rep)
    J <- matrix(1, m.rep, m.rep)
    if (centerfns) {
      V <-
        (1 / (m.rep - 1)) * (t(B) %*% (I_matrix - (1 / m.rep) * J) %*% B)
    } else {
      V <- (1 / (m.rep - 1)) * (t(B) %*% B)
    }
    GCV_score <- 10 ^ 60
    GCVs <- c()
    # Initializes the progress bar
    n_iter <- dim(alpha)[1]
    pb <- txtProgressBar(
      min = 0,
      # Minimum value of the progress bar
      max = n_iter,
      # Maximum value of the progress bar
      style = 3,
      # Progress bar style (also available style = 1 and style = 2)
      width = 50,
      # Progress bar width. Defaults to getOption("width")
      char = "="
    ) # Character used to create the bar
    
    for (j in 1:dim(alpha)[1]) {
      setTxtProgressBar(pb, j)
      I <- I_alpha(mvmfd_obj, alpha[j,])
      D <- I %*% penalty
      L <- t(chol(G + D))
      S <- as.matrix(solve(L))
      rank_cov = qr(S %*% t(G) %*% V %*% G %*% t(S))$rank
      E <- eigen(S %*% t(G) %*% V %*% G %*% t(S))
      if (rank_cov < n) {
        warning(
          "The rank of the coefficient matrix is ",
          rank_cov,
          ". The number of components for computation cannot exceed this rank and has been adjusted accordingly."
        )
        
        n = rank_cov
        u <- E$vectors[, 1:n]
      } else{
        u <- E$vectors
      }
      s_alpha <- sqrtm(solve(G + D))
      s_alpha_tilde <- G_half %*% (solve(G + D)) %*% G_half
      b_temp <- c()
      for (k in 1:n) {
        b_temp <-
          cbind(b_temp, ((t(S) %*% u[, k]) %*% (t(u[, k]) %*% S %*% G %*% t(S) %*% u[, k]) ^
                           (-0.5)))
      }
      v_temp <- B_c %*% G %*% b_temp
      v_temp <- sweep(v_temp, 2, sqrt(diag(t(v_temp) %*% v_temp)), "/")
      GCV_score_temp = gcv_local(
        data = B_c,
        mvmfd_obj = mvmfd_obj,
        G = G,
        G_half = G_half,
        S_smooth = s_alpha_tilde,
        u = v_temp,
        smooth_tuning = alpha[j,]
      )
      GCVs <- c(GCVs, GCV_score_temp)
      
      if (GCV_score_temp < GCV_score) {
        b <- b_temp
        v <- v_temp
        GCV_score <- GCV_score_temp
        GCV_result <- alpha[j,]
      }
    }
    close(pb) # Close the connection
    if (p == 2) {
      GCVs <- matrix(GCVs, nrow = gcv_row, ncol = gcv_column)
    }
    temp_count <- 0
    pc <- list()
    for (i in 1:p) {
      index_start <- (temp_count + 1)
      index_end <- (temp_count + prod(mvmfd_obj$basis$nbasis[[i]]))
      pc[[i]] <- b[index_start:index_end,]
      temp_count <- temp_count + prod(mvmfd_obj$basis$nbasis[[i]])
    }
    variance <- diag(t(b) %*% G %*% V %*% G %*% b)
    lsv <-
      (B_c %*% G %*% b) %*% solve(diag(sqrt(diag(
        t(B_c %*% G %*% b) %*% (B_c %*% G %*% b)
      )), ncol = ncol(b)))
    bbbb <- c()
    for (k in 1:p) {
      bbbb <- rbind(bbbb, pc[[k]])
    }
    return(list(pc, lsv, variance, GCV_result, GCVs))
  }
