#' @importFrom expm sqrtm
#' @importFrom utils  txtProgressBar setTxtProgressBar
#'
#sparse penalty function
sparse_pen_fun <- function(y, tuning_parameter, type, alpha = 3.7) {
  y_sorted <- sort(abs(y))
  lambda = y_sorted[tuning_parameter]
  if (tuning_parameter == 0) {
    return(y)
  }
  if (type == "soft") {
    return(sign(y) * pmax(abs(y) - lambda, 0))
  }
  else if (type == "hard") {
    return(ifelse(abs(y) > lambda, y, 0))
  }
  else if (type == "SCAD") {
    res <- ifelse(
      abs(y) <= 2 * lambda,
      sign(y) * pmax(abs(y) - lambda, 0),
      ifelse(
        abs(y) <= alpha * lambda,
        ((alpha - 1) * y - sign(y) * alpha * lambda) / (alpha - 2),
        y
      )
    )
    return(res)
  }
}

#sequential power algorithm
init_sequential = function(data,
                           sparse_tuning_result,
                           sparse_tuning_type,
                           S_smooth = NULL,
                           S_2_inverse = NULL,
                           G_half_inverse = NULL,
                           G_half = NULL,
                           cv_flag = FALSE) {
  v_old = svd(as.matrix(data))$v[, 1]
  errors = 10 ^ 60
  while (errors > 0.00001) {
    u_new = sparse_pen_fun(y = data %*% v_old,
                           tuning_parameter = sparse_tuning_result,
                           sparse_tuning_type)
    u_new = u_new / norm(u_new, type = "2")
    if (cv_flag == TRUE) {
      #incorporate power algorithm in CV sparse tuning selection
      v_new = t(data) %*% u_new
      v_new = v_new / norm(v_new, type = "2")
      errors = sum((v_new - v_old) ^ 2)
      v_old = v_new
    } else{
      v_new = S_smooth %*% t(data) %*% u_new
      errors = sum((v_new - v_old) ^ 2)
      v_old = v_new
    }
  }
  if (cv_flag == TRUE) {
    return(u_new)
  }
  else{
    v_new = G_half_inverse %*% v_new
    v_new = v_new %*% solve(sqrt(t(v_new) %*% S_2_inverse %*% v_new))
    return(list(v_new, u_new))
  }
}

#joint power for smoothing
init_joint = function(data,
                      S_smooth = NULL,
                      S_2_inverse = NULL,
                      G_half_inverse = NULL,
                      G_half = NULL,
                      n = n) {
  v_old = svd(data)$v[, 1:n]
  errors = 10 ^ 60
  while (errors > 0.00001) {
    u_new = data %*% v_old
    u_new = sweep(u_new, 2, sqrt(diag(t(u_new) %*% u_new)), "/")
    v_new = S_smooth %*% qr.Q(qr(as.matrix(t(data) %*% u_new)))
    errors = sum((v_new - v_old) ^ 2) / n
    v_old = v_new
  }
  v_new = G_half_inverse %*% v_new
  v_new = sweep(v_new, 2, sqrt(diag(t(v_new) %*% S_2_inverse %*% v_new)), "/")
  return(list(v_new, u_new))
}


#computing cv score for sparse tuning
cv_local = function(data,
                    G_half,
                    K_fold,
                    sparse_tuning_single,
                    sparse_tuning_type,
                    shuffled_row,
                    group_size) {
  data_double_tilde = t(data %*% G_half)
  error_score_sparse = 0
  for (k in 1:K_fold) {
    rows_to_remove = shuffled_row[((k - 1) * group_size + 1):((k) * group_size)]
    data_train = data_double_tilde[-rows_to_remove,]
    data_test = data_double_tilde[rows_to_remove,]
    u_test = init_sequential(t(data_train),
                             sparse_tuning_single,
                             sparse_tuning_type,
                             cv_flag = TRUE)
    v_test = data_test %*% u_test
    v_test_smooth_back = (data_double_tilde %*% u_test)[rows_to_remove,]
    data_test_smooth_back = t(data_double_tilde)[, rows_to_remove]
    error_score_sparse = error_score_sparse + sum((
      t(data_test_smooth_back) - v_test_smooth_back %*% t(u_test)
    ) ^ 2)
  }
  return(error_score_sparse / ncol(data))
}

#computing gcv score for smoothing tuning
gcv_local = function(data,
                     mvmfd_obj,
                     G,
                     G_half,
                     S_smooth,
                     u,
                     smooth_tuning) {
  p = mvmfd_obj$nvar
  indices <-
    sapply(1:p, function(i)
      prod(mvmfd_obj$basis$nbasis[[i]]))
  C_subtilde = data %*% G_half
  if (all(smooth_tuning == 0)) {
    error_smooth_score <- 0
  } else {
    error_smooth_score <- 0
    start_index <- 1
    
    for (i in 1:p) {
      end_index <- start_index + indices[i] - 1
      
      if (smooth_tuning[i] == 0) {
        error_smooth_score_i <- 0
      } else {
        s_alpha_tilde_i <-
          S_smooth[start_index:end_index, start_index:end_index]
        C_subtilde_i <- C_subtilde[, start_index:end_index]
        error_smooth_score_i <-
          sum(((diag(indices[i]) - s_alpha_tilde_i) %*% (t(
            C_subtilde_i
          ) %*% u)) ^ 2) /
          (1 - sum(diag(s_alpha_tilde_i)) / indices[i]) ^ 2
      }
      
      error_smooth_score <-
        error_smooth_score + error_smooth_score_i
      start_index <- end_index + 1
    }
  }
  
  return(error_smooth_score)
}

# Function to handle smooth tuning selection with progress bar and GCV score calculation
handle_smooth_tuning <-
  function(data,
           G_half,
           G,
           S_smooth,
           S_2_inverse,
           G_half_inverse,
           mvmfd_obj,
           sparse_tuning_selection = NULL,
           sparse_tuning_type = NULL,
           smooth_tuning,
           CV_score_smooth,
           power_type = "sequential",
           n = NULL,
           pb,
           count) {
    gcv_scores <- NULL
    smooth_tuning_selection <- NULL
    index_selection <- NULL
    if (is.null(smooth_tuning)) {
      count <- count + 1
      setTxtProgressBar(pb, count)
      smooth_tuning_selection <-
        expand.grid(lapply(rep(0, mvmfd_obj$nvar), function(x)
          x[1]))
      index_selection <- 1
    } else {
      for (smooth_index in 1:dim(smooth_tuning)[1]) {
        count <- count + 1
        setTxtProgressBar(pb, count)
        if (all(smooth_tuning == 0)) {
          S_smooth[[smooth_index]] <- diag(dim(G)[1])
        }
        if (power_type == "sequential") {
          test_temp <-
            init_sequential(
              data %*% G_half,
              sparse_tuning_selection,
              sparse_tuning_type,
              S_smooth[[smooth_index]],
              S_2_inverse[[smooth_index]],
              G_half_inverse,
              G_half
            )
        } else {
          test_temp <-
            init_joint(data %*% G_half,
                       S_smooth[[smooth_index]],
                       S_2_inverse[[smooth_index]],
                       G_half_inverse,
                       G_half,
                       n = n)
        }
        v_temp <- test_temp[[2]]
        smooth_score <-
          gcv_local(data,
                    mvmfd_obj,
                    G,
                    G_half,
                    S_smooth[[smooth_index]],
                    v_temp,
                    smooth_tuning = smooth_tuning[smooth_index,])
        gcv_scores <- c(gcv_scores, smooth_score)
        if (smooth_score <= CV_score_smooth) {
          CV_score_smooth <- smooth_score
          smooth_tuning_selection <- smooth_tuning[smooth_index,]
          index_selection <- smooth_index
        }
      }
    }
    return(
      list(
        smooth_tuning_selection = smooth_tuning_selection,
        index_selection = index_selection,
        gcv_scores = gcv_scores
      )
    )
  }

# Function to handle sparse tuning selection
handle_sparse_tuning <-
  function(data,
           G_half,
           sparse_tuning,
           sparse_tuning_type,
           K_fold,
           shuffled_row,
           group_size,
           CV_score_sparse,
           pb) {
    count <- 0
    cv_scores <- c()
    sparse_tuning_selection <- NULL
    if (is.null(sparse_tuning)) {
      count <- count + 1
      setTxtProgressBar(pb, count)
      cv_scores <- NULL
      sparse_tuning_selection <- 0
    } else {
      for (sparse_tuning_single in sparse_tuning) {
        count <- count + 1
        setTxtProgressBar(pb, count)
        sparse_score <-
          cv_local(
            data,
            G_half,
            K_fold,
            sparse_tuning_single,
            sparse_tuning_type,
            shuffled_row,
            group_size
          )
        cv_scores <- c(cv_scores, sparse_score)
        if (sparse_score <= CV_score_sparse) {
          CV_score_sparse <- sparse_score
          sparse_tuning_selection <- sparse_tuning_single
        }
      }
    }
    return(
      list(
        sparse_tuning_selection = sparse_tuning_selection,
        cv_scores = cv_scores,
        CV_score_sparse = CV_score_sparse
      )
    )
  }

# Function for cv_gcv_sequential
cv_gcv_sequential <-
  function(data,
           mvmfd_obj,
           smooth_tuning,
           sparse_tuning,
           sparse_tuning_type,
           K_fold,
           G,
           G_half,
           G_half_inverse,
           S_smooth,
           S_2_inverse) {
    CV_score_sparse <- CV_score_smooth <- Inf
    result <- c()
    count <- 0
    shuffled_row <- sample(ncol(data))
    group_size <- length(shuffled_row) / K_fold
    n_iter <-
      (if (is.null(smooth_tuning))
        1
       else
         dim(smooth_tuning)[1]) + (if (is.null(sparse_tuning))
           1
           else
             length(sparse_tuning))
    pb <-
      txtProgressBar(
        min = 0,
        max = n_iter,
        style = 3,
        width = 50,
        char = "="
      )
    # Handle sparse tuning
    sparse_tuning_result <-
      handle_sparse_tuning(
        data,
        G_half,
        sparse_tuning,
        sparse_tuning_type,
        K_fold,
        shuffled_row,
        group_size,
        CV_score_sparse,
        pb = pb
      )
    sparse_tuning_selection <-
      sparse_tuning_result$sparse_tuning_selection
    cv_scores <- sparse_tuning_result$cv_scores
    CV_score_sparse <- sparse_tuning_result$CV_score_sparse
    # Handle smooth tuning
    smooth_tuning_result <-
      handle_smooth_tuning(
        data,
        G_half,
        G,
        S_smooth,
        S_2_inverse,
        G_half_inverse,
        mvmfd_obj,
        sparse_tuning_selection,
        sparse_tuning_type,
        smooth_tuning,
        CV_score_smooth,
        power_type = "sequential",
        pb = pb,
        count = (if (is.null(sparse_tuning))
          1
          else
            length(sparse_tuning))
      )
    smooth_tuning_selection <-
      smooth_tuning_result$smooth_tuning_selection
    index_selection <- smooth_tuning_result$index_selection
    gcv_scores <- smooth_tuning_result$gcv_scores
    
    close(pb)
    result <-
      list(
        sparse_tuning_selection,
        smooth_tuning_selection,
        index_selection,
        cv_scores,
        gcv_scores
      )
    return(result)
  }

# Function for gcv_joint
gcv_joint <-
  function(data,
           G_half,
           G,
           S_smooth,
           S_2_inverse,
           G_half_inverse,
           mvmfd_obj,
           smooth_tuning,
           n) {
    CV_score_smooth <- Inf
    pb <-
      txtProgressBar(
        min = 0,
        max = dim(smooth_tuning)[1],
        style = 3,
        width = 50,
        char = "="
      )
    smooth_tuning_result <-
      handle_smooth_tuning(
        data,
        G_half,
        G,
        S_smooth,
        S_2_inverse,
        G_half_inverse,
        mvmfd_obj,
        smooth_tuning = smooth_tuning,
        CV_score_smooth = CV_score_smooth,
        power_type = "joint",
        n = n,
        pb = pb,
        count = 0
      )
    smooth_tuning_selection <-
      smooth_tuning_result$smooth_tuning_selection
    index_selection <- smooth_tuning_result$index_selection
    gcv_scores <- smooth_tuning_result$gcv_scores
    
    close(pb)
    return(
      list(
        smooth_tuning_selection = smooth_tuning_selection,
        index_selection = index_selection,
        gcv_scores = gcv_scores
      )
    )
  }


ordinal_msg <- function(i) {
  if (i == 1) {
    return(paste0(i, "st"))
  } else if (i == 2) {
    return(paste0(i, "nd"))
  } else if (i == 3) {
    return(paste0(i, "rd"))
  }
  else {
    return(paste0(i, "th"))
  }
}

# Define a function to process mvmfd_obj
centralized_mvmfd <- function(mvmfd_obj, centerfns = TRUE) {
  p <- mvmfd_obj$nvar
  C <- c()
  
  for (i in 1:p) {
    coefs_i <- mvmfd_obj$coefs[[i]]
    
    if (is.matrix(coefs_i)) {
      c <- if (centerfns)
        coefs_i - rowMeans(coefs_i)
      else
        coefs_i
    } else {
      cc <- apply(coefs_i, 3, as.vector)
      c <- if (centerfns)
        cc - rowMeans(cc)
      else
        cc
    }
    
    C <- rbind(C, c)
  }
  
  return(C)
}

# Function to handle variance calculation and update
handle_variance_update <-
  function(i,
           n,
           C,
           G,
           v_total,
           mvmfd_obj,
           all_equal_check,
           sparse_tuning,
           pc,
           lsv,
           v,
           u,
           G_half,
           test_result,
           temp_count,
           p) {
    for (j in 1:p) {
      index_start <- (temp_count + 1)
      index_end <- (temp_count + prod(mvmfd_obj$basis$nbasis[[j]]))
      if (i == 1) {
        pc[[j]] <- v[index_start:index_end,]
      } else {
        pc[[j]] <- cbind(pc[[j]], v[index_start:index_end,])
      }
      temp_count <- index_end
    }
    lsv = cbind(lsv, u)
    v_total = cbind(v_total, v)
    
    if (i == 1 ||
        all(all_equal_check) &
        (is.null(sparse_tuning) || all(unique(sparse_tuning) == 0))) {
      CGv <- C %*% G %*% v
      variance <- t(CGv) %*% CGv / (mvmfd_obj$nobs - 1)
    } else {
      CGv <- C %*% G %*% v_total
      G_pc = t(v_total) %*% G %*% v_total
      coef_pc = CGv %*% solve(G_pc)
      total_variance = sum(diag((coef_pc %*% t(v_total)) %*% G %*% t(coef_pc %*%
                                                                       t(v_total))))
      G_pc_pre = t(v_total[,-i]) %*% G %*% v_total[,-i]
      coef_pc_pre = CGv[,-i] %*% solve(G_pc_pre)
      total_variance_previous = sum(diag((coef_pc_pre %*% t(v_total[,-i])) %*% G %*% t(coef_pc_pre %*%
                                                                                         t(v_total[,-i]))))
      variance = (total_variance - total_variance_previous) / (mvmfd_obj$nobs - 1)
    }
    return(list(
      pc = pc,
      lsv = lsv,
      v_total = v_total,
      variance = variance
    ))
  }

sequential_power <-
  function(mvmfd_obj,
           n,
           smooth_tuning,
           smooth_tuning_type,
           sparse_tuning,
           sparse_tuning_type,
           centerfns,
           alpha_orth,
           K_fold,
           sparse_CV,
           smooth_GCV) {
    p <- mvmfd_obj$nvar
    smooth_penalty <-
      pen_fun(mvmfd_obj, type = smooth_tuning_type)
    
    #######centralize########
    C <- centralized_mvmfd(mvmfd_obj, centerfns)
    ########some initial setting#######
    C <- t(C)
    lsv <- c()
    pc <- list()
    variance <- vector()
    smooth_tuning_result  <- list()
    sparse_tuning_result  <- list()
    G <- as.matrix(mvmfd_obj$basis$gram)
    G_half <- sqrtm(G)
    G_half_inverse = solve(G_half)
    all_equal_check <-
      sapply(smooth_tuning, function(x)
        length(unique(x)) == 1)
    rank_C = qr(C %*% G_half)$rank
    if (rank_C < n) {
      warning(
        "The rank of the coefficient matrix is ",
        rank_C,
        ". The number of components for computation cannot exceed this rank and has been adjusted accordingly."
      )
      n = rank_C
    }
    #########matrix input of smoothing parameters###########
    if (smooth_GCV == FALSE) {
      v_total = c()
      # I_a <- D <- S_2 <- S_smooth <- S_2_inverse <- list()
      S_smooth <- S_2_inverse <- list()
      GCV_score = c()
      if (sparse_CV == FALSE) {
        CV_score = c()
      } else{
        CV_score = list()
      }
      for (i in 1:n) {
        cat(sprintf("Computing the %s PC...\n", ordinal_msg(i)))
        if (is.null(smooth_tuning)) {
          smooth_tuning_temp = expand.grid(lapply(rep(0, mvmfd_obj$nvar), function(x)
            x[1]))
        } else{
          smooth_tuning_temp = expand.grid(lapply(smooth_tuning, function(x)
            x[i]))
        }
        if (i == 1) {
          C_temp = C
        }
        else{
          b_original = t(C_temp) %*% u
          C_temp = C_temp - u %*% t(b_original)
        }
        D <-
          I_alpha(mvmfd_obj, smooth_tuning_temp) %*% smooth_penalty
        S_2 <- solve(G + D)
        S_2_inverse[[1]] = solve(S_2)
        S_smooth[[1]] <- G_half %*% (S_2) %*% G_half
        
        if (!is.null(sparse_tuning)) {
          sparse_tuning_temp <-
            if (sparse_CV == FALSE)
              sparse_tuning[i]
          else
            sparse_tuning
        }
        cv_result = cv_gcv_sequential(
          data = C_temp,
          mvmfd_obj = mvmfd_obj,
          smooth_tuning = if (is.null(smooth_tuning))
            smooth_tuning
          else
            smooth_tuning_temp,
          sparse_tuning = if (is.null(sparse_tuning))
            sparse_tuning
          else
            sparse_tuning_temp,
          sparse_tuning_type = sparse_tuning_type,
          K_fold = K_fold,
          G = G,
          G_half = G_half,
          G_half_inverse = G_half_inverse,
          S_smooth = S_smooth,
          S_2_inverse = S_2_inverse
        )
        sparse_result = cv_result[[1]]
        smooth_result_index = cv_result[[3]]
        if (sparse_CV == FALSE) {
          CV_score = c(CV_score, cv_result[[4]])
        } else{
          CV_score[[i]] = cv_result[[4]]
        }
        GCV_score = c(GCV_score, cv_result[[5]])
        test_result = init_sequential(
          C_temp %*% G_half,
          sparse_result,
          sparse_tuning_type,
          S_smooth[[1]],
          S_2_inverse[[1]],
          G_half_inverse,
          G_half
        )
        
        u = test_result[[2]]
        v = test_result[[1]]
        smooth_tuning_result[[i]] = smooth_tuning_temp[smooth_result_index,]
        sparse_tuning_result[[i]] = sparse_result
        temp_count <- 0
        variance_result <-
          handle_variance_update(
            i,
            n,
            C,
            G,
            v_total,
            mvmfd_obj,
            all_equal_check,
            sparse_tuning,
            pc,
            lsv,
            v,
            u,
            G_half,
            test_result,
            temp_count,
            p
          )
        pc <- variance_result$pc
        lsv <- variance_result$lsv
        v_total <- variance_result$v_total
        variance[i] <- variance_result$variance
      }
    }
    #########sequential inputs of smoothing parameters###########
    else{
      if (is.null(smooth_tuning)) {
        smooth_tuning_temp = expand.grid(lapply(rep(0, mvmfd_obj$nvar), function(x)
          x[1]))
      } else{
        smooth_tuning_temp <- expand.grid(smooth_tuning)
      }
      
      S_smooth <- S_2_inverse <- list()
      cat("Preprocessing...\n")
      n_iter1 <- dim(smooth_tuning_temp)[1]
      pb <-
        txtProgressBar(
          min = 0,
          # Minimum value of the progress bar
          max = n_iter1,
          # Maximum value of the progress bar
          style = 3,
          # Progress bar style (also available style = 1 and style = 2)
          width = 50,
          # Progress bar width. Defaults to getOption("width")
          char = "="
        )   # Character used to create the bar
      #####Do the following computation in advance to save computational cost#####
      for (smooth_index in 1:dim(smooth_tuning_temp)[1]) {
        setTxtProgressBar(pb, smooth_index)
        D <-
          I_alpha(mvmfd_obj, smooth_tuning_temp[smooth_index,]) %*% smooth_penalty
        S_2 <- solve(G + D)
        S_2_inverse[[smooth_index]] = solve(S_2)
        S_smooth[[smooth_index]] <- G_half %*% (S_2) %*% G_half
      }
      close(pb)
      v_total = c()
      GCV_score = list()
      if (sparse_CV == FALSE) {
        CV_score = c()
      } else{
        CV_score = list()
      }
      for (i in 1:n) {
        cat(sprintf("Computing the %s PC...\n", ordinal_msg(i)))
        if (i == 1) {
          C_temp = C
        }
        else{
          b_original = t(C_temp) %*% u
          C_temp = C_temp - u %*% t(b_original)
        }
        if (!is.null(sparse_tuning)) {
          sparse_tuning_temp <-
            if (sparse_CV == FALSE)
              sparse_tuning[i]
          else
            sparse_tuning
        }
        cv_result = cv_gcv_sequential(
          data = C_temp,
          mvmfd_obj = mvmfd_obj,
          smooth_tuning = if (is.null(smooth_tuning))
            smooth_tuning
          else
            smooth_tuning_temp,
          sparse_tuning = if (is.null(sparse_tuning))
            sparse_tuning
          else
            sparse_tuning_temp,
          sparse_tuning_type = sparse_tuning_type,
          K_fold = K_fold,
          G = G,
          G_half = G_half,
          G_half_inverse = G_half_inverse,
          S_smooth = S_smooth,
          S_2_inverse = S_2_inverse
        )
        sparse_result = cv_result[[1]]
        smooth_result_index = cv_result[[3]]
        if (sparse_CV == FALSE) {
          CV_score = c(CV_score, cv_result[[4]])
        } else{
          CV_score[[i]] = cv_result[[4]]
        }
        GCV_score[[i]] = cv_result[[5]]
        test_result = init_sequential(
          C_temp %*% G_half,
          sparse_result,
          sparse_tuning_type,
          S_smooth[[smooth_result_index]],
          S_2_inverse[[smooth_result_index]],
          G_half_inverse,
          G_half
        )
        u = test_result[[2]]
        v = test_result[[1]]
        smooth_tuning_result[[i]] = smooth_tuning_temp[smooth_result_index,]
        sparse_tuning_result[[i]] = sparse_result
        temp_count <- 0
        variance_result <-
          handle_variance_update(
            i,
            n,
            C,
            G,
            v_total,
            mvmfd_obj,
            all_equal_check,
            sparse_tuning,
            pc,
            lsv,
            v,
            u,
            G_half,
            test_result,
            temp_count,
            p
          )
        pc <- variance_result$pc
        lsv <- variance_result$lsv
        v_total <- variance_result$v_total
        variance[i] <- variance_result$variance
      }
      if (is.null(smooth_tuning)) {
        GCV_score = NULL
      }
      if (is.null(sparse_tuning)) {
        CV_score = NULL
      }
    }
    return(
      list(
        pc,
        lsv,
        variance,
        smooth_tuning_result,
        sparse_tuning_result,
        CV_score,
        GCV_score
      )
    )
  }


# joint smooth and sparse power algorithm
joint_power <-
  function(mvmfd_obj,
           n,
           smooth_tuning,
           smooth_tuning_type,
           centerfns,
           alpha_orth) {
    p <- mvmfd_obj$nvar
    smooth_penalty <-
      pen_fun(mvmfd_obj, type = smooth_tuning_type)
    
    #######centralize########
    C <- centralized_mvmfd(mvmfd_obj, centerfns)
    #########################
    
    ########some initial setting#######
    C <- t(C)
    lsv <- c()
    pc <- list()
    variance <- vector()
    smooth_tuning_result  <- list()
    sparse_tuning_result  <- list()
    G <- as.matrix(mvmfd_obj$basis$gram)
    G_half <- sqrtm(G)
    G_half_inverse = solve(G_half)
    ###################################
    
    rank_C = qr(C %*% G_half)$rank
    if (rank_C < n) {
      warning(
        "The rank of the coefficient matrix is ",
        rank_C,
        ". The number of components for computation cannot exceed this rank and has been adjusted accordingly."
      )
      n = rank_C
    }
    
    #####smoothing parameter#######
    if (is.null(smooth_tuning)) {
      smooth_tuning_temp = expand.grid(lapply(rep(0, mvmfd_obj$nvar), function(x)
        x[1]))
    } else{
      smooth_tuning_temp <- expand.grid(smooth_tuning)
    }
    ####################################
    
    ####joint power####
    # I_a <- D <- S_2 <- S_smooth <- S_2_inverse <- list()
    S_smooth <- S_2_inverse <- list()
    cat("Preprocessing...\n")
    n_iter1 <- dim(smooth_tuning_temp)[1]
    pb <-
      txtProgressBar(
        min = 0,
        # Minimum value of the progress bar
        max = n_iter1,
        # Maximum value of the progress bar
        style = 3,
        # Progress bar style (also available style = 1 and style = 2)
        width = 50,
        # Progress bar width. Defaults to getOption("width")
        char = "="
      )   # Character used to create the bar
    #####Do the following computation in advance to save computational cost#####
    for (smooth_index in 1:dim(smooth_tuning_temp)[1]) {
      setTxtProgressBar(pb, smooth_index)
      D <-
        I_alpha(mvmfd_obj, smooth_tuning_temp[smooth_index,]) %*% smooth_penalty
      S_2 <- solve(G + D)
      S_2_inverse[[smooth_index]] = solve(S_2)
      S_smooth[[smooth_index]] <- G_half %*% (S_2) %*% G_half
    }
    close(pb)
    cat(sprintf("Computing PCs...\n"))
    # cv_result = gcv_joint(data = C, mvmfd_obj = mvmfd_obj, smooth_tuning = smooth_tuning, G = G, G_half = G_half, G_half_inverse = G_half_inverse, S_smooth = S_smooth, S_2_inverse = S_2_inverse, n = n)
    cv_result = gcv_joint(
      data = C,
      mvmfd_obj = mvmfd_obj,
      smooth_tuning = if (is.null(smooth_tuning))
        smooth_tuning
      else
        smooth_tuning_temp,
      G = G,
      G_half = G_half,
      G_half_inverse = G_half_inverse,
      S_smooth = S_smooth,
      S_2_inverse = S_2_inverse,
      n = n
    )
    smooth_result_index = cv_result[[2]]
    GCV_score = cv_result[[3]]
    test_result = init_joint(C %*% G_half,
                             S_smooth[[smooth_result_index]],
                             S_2_inverse[[smooth_result_index]],
                             G_half_inverse,
                             G_half,
                             n = n)
    u = test_result[[2]]
    v = test_result[[1]]
    smooth_tuning_result = smooth_tuning_temp[smooth_result_index,]
    temp_count <- 0
    
    temp_count <- 0
    for (j in 1:p) {
      index_start <- (temp_count + 1)
      index_end <- (temp_count + prod(mvmfd_obj$basis$nbasis[[j]]))
      pc[[j]] <- v[index_start:index_end,]
      temp_count <- index_end
    }
    
    lsv = cbind(lsv, u)
    v_total = v
    for (k in 1:n) {
      CGv <- C %*% G %*% v[, k]
      variance[k] <- t(CGv) %*% CGv / (mvmfd_obj$nobs - 1)
    }
    return(list(pc, lsv, variance, smooth_tuning_result, GCV_score))
  } 