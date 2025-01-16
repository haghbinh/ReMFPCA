#' @title A Class for `ReMFPCA` objects
#'
#' @description
#' The `remfpca` class represents regularized functional principal components components.
#'

#' @field pc_mfd An object of class `mvmfd` where the first indices (fields)
#' represents harmonics and  second indices represents variables
#' @field lsv = Left singular values vectors
#' @field values = The set of eigenvalues
#' @field smooth_tuning = The list of smoothing penalties parameters
#' @field sparse_tuning = The list of sparse penalties parameters
#' @field GCVs = Generalized cross validations scores of smoothing penalties parameters.
#'               If both smoothing and sparse tuning penalties are used in the ReMFPCA method,
#'               this represents the conditional generalized cross-validation scores, which
#'               means it is computed based on the optimal sparse tuning parameter selected via cross validation.
#' @field CVs = Cross validations scores of sparse penalties parameters
#' @field mean_mfd A multivariate functional data object giving the mean function

#' @examples
#' require(fda)
#' # Brownian Bridge simulation on [0,1]
#' M <- 110 # number of components
#' N <- 20 # number of instances
#' n <- 100 # number of grides
#' t0 <- seq(0, 1, len = n)
#' j <- 1:M
#' alpha1 <- list(a1 = 2^seq(0, 1, length.out = 3), a2 = 2^seq(0, 1, length.out = 3))
#' sparse_tuning = as.integer(seq(1, N-1, length.out = 10))
#' psi_1 <- function(t, m) sin(m * pi * t) # eigenfunction of BB
#' psi_2 <- function(t, m) sin((2 * m - 1) * pi / 2 * t) # eigenfunction of BM
#' PC_1 <- outer(t0, j, FUN = psi_1) # n by M matrix
#' PC_2 <- outer(t0, j, FUN = psi_2) # n by M matrix
#' Z <- matrix(rnorm(N * M), nr = M)
#' lambda <- matrix(2 / (pi * (2 * j - 1)), nr = M, nc = N)
#' X_1t <- PC_1 %*% (lambda * Z)
#' X_2t <- PC_2 %*% (lambda * Z)
#' noise <- rnorm(n * N, 0, 0.1)
#' X_1 <- X_1t + noise
#' X_2 <- X_2t + noise
#' bs <- create.bspline.basis(c(0, 1), 51)
#' mdbs <- Basismfd(bs)
#' mfd1 <- Mfd(X = X_1, mdbs = mdbs)
#' mfd2 <- Mfd(X = X_2, mdbs = mdbs)
#' mvmfd_obj <- Mvmfd(mfd1, mfd2)
#' k <- 2
#' # Non Regularized MFPCA based on sequential power algorithm
#' Re0 <- Remfpca(mvmfd_obj, ncomp = k, smooth_GCV = FALSE, sparse_CV = FALSE)
#' fpc0 <- Re0$pc_mfd
#' scores0 <- inprod_mvmfd(mvmfd_obj, fpc0)
#' dim(scores0)
#' # Smooth MFPCA based on sequential power algorithm
#' Re1 <- Remfpca(mvmfd_obj, ncomp = k, smooth_tuning = alpha1)
#' # Smooth and sparse MFPCA based on sequential power algorithm
#' Re2 <- Remfpca(mvmfd_obj, ncomp = k, smooth_tuning = alpha1, sparse_tuning = sparse_tuning)
#' # Smooth MFPCA based on joint power algorithm
#' Re3 <- Remfpca(mvmfd_obj, ncomp = k, smooth_tuning = alpha1, alpha_orth = TRUE)
#' # Smooth MFPCA based on eigen decomposition algorithm
#' Re4 <- Remfpca(mvmfd_obj, ncomp = k, smooth_tuning = alpha1, method = "eigen")


#' @import R6
#' @importFrom fda getbasispenalty
#' @importFrom Matrix bdiag
#' @importFrom expm sqrtm
#' @export
remfpca <- R6::R6Class(
  "remfpca",
  public = list(
    #' @description
    #' Initialize the `remfpca` class.
    #' @param mvmfd_obj An `mvmfd` object representing the multivariate functional data.
    #' @param method A character string specifying the approach to be used for MFPCA computation.
    #'               Options are "power" (the default) or "eigen".
    #' @param ncomp The number of functional principal components to retain.
    #' @param smooth_tuning A list or vector specifying the smoothing regularization parameter(s).
    #' @param sparse_tuning A list or vector specifying the sparsity regularization parameter(s).
    #' @param centerfns Logical. Whether to center the functional data before analysis.
    #' @param alpha_orth Logical. Whether to perform orthogonalization of the regularization parameters.
    #' @param smoothing_type The type of smoothing penalty to be applied.
    #' @param sparse_type The type of sparse penalty to be applied.
    #' @param K_fold An integer specifying the number of folds in cross-validation.
    #' @param sparse_CV Logical. Whether cross-validation should be applied for sparse tuning.
    #' @param smooth_GCV Logical. Whether generalized cross-validation should be applied for smoothing tuning.
    
    # initialize = function(mvmfd_obj, method = "power", ncomp, smooth_tuning = NULL, sparse_tuning = 0, centerfns = TRUE, alpha_orth = FALSE, smoothing_type = "coefpen", sparse_type = "soft", K_fold = 30, sparsity_CV = "marginal") {
    initialize = function(mvmfd_obj,
                          method = "power",
                          ncomp,
                          smooth_tuning = NULL,
                          sparse_tuning = NULL,
                          centerfns = TRUE,
                          alpha_orth = FALSE,
                          smoothing_type = "coefpen",
                          sparse_type = "soft",
                          K_fold = 30,
                          sparse_CV,
                          smooth_GCV) {
      # if (is.numeric(smooth_tuning)) smooth_tuning <- as.list(smooth_tuning)
      # if (is.vector(smooth_tuning)) smooth_tuning <- as.list(smooth_tuning)
      # if (is.vector(smooth_tuning)&& !is.list(smooth_tuning)) smooth_tuning <- list(smooth_tuning)
      if (is.mfd(mvmfd_obj))
        mvmfd_obj <- mvmfd$new(mvmfd_obj)
      
      # if (ncomp > Reduce(`+`, mvmfd_obj$basis$nbasis)) {
      #   ncomp <- Reduce(`+`, mvmfd_obj$basis$nbasis)
      #   warning("The maximum number of components for computation should not exceed the total number of basis functions across all variables. The number of components has been set to the total number of basis functions.")
      # }
      
      if (method == "power" &
          alpha_orth == "FALSE") {
        # Adjust the vector length to match the required dimensions if they are incorrect
        if (is.vector(smooth_tuning) &
            !is.list(smooth_tuning)) {
          if (smooth_GCV == FALSE) {
            if (length(smooth_tuning) != ncomp) {
              warning(
                "The length of 'smooth_tuning' did not match 'ncomp' and has been adjusted accordingly.",
                call. = FALSE
              )
              smooth_tuning <-
                rep(smooth_tuning, length.out = ncomp)
            }
            smooth_tuning <-
              replicate(mvmfd_obj$nvar, smooth_tuning, simplify = FALSE)
          }
          else{
            warning(
              "The length of 'smooth_tuning' did not match 'mvmfd_obj$nvar' and has been adjusted accordingly.",
              call. = FALSE
            )
            smooth_tuning <-
              replicate(mvmfd_obj$nvar, smooth_tuning, simplify = FALSE)
          }
        }
        
        # Adjust the matrix to match the required dimensions if they are incorrect
        else if (is.matrix(smooth_tuning)) {
          if (smooth_GCV == FALSE) {
            if (!all(dim(smooth_tuning) == c(mvmfd_obj$nvar, ncomp))) {
              smooth_tuning <-
                smooth_tuning[rep(1:nrow(smooth_tuning), length.out = mvmfd_obj$nvar), rep(1:ncol(smooth_tuning), length.out = ncomp)]
              # print(smooth_tuning)
              smooth_tuning <-
                split(smooth_tuning, row(smooth_tuning))
              # print(smooth_tuning)
              warning(
                "The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.",
                call. = FALSE
              )
            } else{
              smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
            }
          }
          else{
            if (dim(smooth_tuning)[1] != mvmfd_obj$nvar) {
              smooth_tuning <-
                smooth_tuning[rep(1:nrow(smooth_tuning), length.out = mvmfd_obj$nvar), , drop = FALSE][1:mvmfd_obj$nvar, , drop = FALSE]
              smooth_tuning <-
                split(smooth_tuning, row(smooth_tuning))
              warning(
                "The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.",
                call. = FALSE
              )
            }
            else{
              smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
            }
          }
        }
        
        # Adjust the list length and element sizes to match the required dimensions if they are incorrect
        else if (is.list(smooth_tuning)) {
          if (smooth_GCV == FALSE) {
            if (length(smooth_tuning) != mvmfd_obj$nvar) {
              warning("Adjusting 'smooth_tuning' to match 'mvmfd_obj$nvar'.",
                      call. = FALSE)
              smooth_tuning <-
                rep(smooth_tuning, length.out = mvmfd_obj$nvar)
            }
            smooth_tuning <-
              lapply(smooth_tuning, function(vec) {
                if (length(vec) != ncomp) {
                  warning("Adjusting vector length in 'smooth_tuning' to match 'ncomp'.",
                          call. = FALSE)
                  vec <-
                    rep(vec, length.out = ncomp)
                }
                vec
              })
          }
          else{
            if (length(smooth_tuning) != mvmfd_obj$nvar) {
              warning("Adjusting 'smooth_tuning' to match 'mvmfd_obj$nvar'.",
                      call. = FALSE)
              smooth_tuning <-
                rep(smooth_tuning, length.out = mvmfd_obj$nvar)
            }
            
            # smooth_tuning <- rep(smooth_tuning, length.out = mvmfd_obj$nvar)[1:mvmfd_obj$nvar]
          }
        }
        
        if (!is.null(smooth_tuning)) {
          names(smooth_tuning) <- paste0("var", 1:mvmfd_obj$nvar)
        }
        
        # Adjust the list length and element sizes to match the required dimensions if they are incorrect
        if (sparse_CV == FALSE &
            length(sparse_tuning) != ncomp & !is.null(sparse_tuning)) {
          warning(
            "The length of 'sparse_tuning' did not match 'ncomp' and has been adjusted accordingly.",
            call. = FALSE
          )
          sparse_tuning <-
            rep(sparse_tuning, length.out = ncomp)
        }
        
        result <-
          sequential_power(
            mvmfd_obj = mvmfd_obj,
            n = ncomp,
            smooth_tuning = smooth_tuning,
            sparse_tuning = sparse_tuning,
            centerfns = centerfns,
            alpha_orth = alpha_orth,
            smooth_tuning_type = smoothing_type,
            sparse_tuning_type = sparse_type,
            K_fold = K_fold,
            sparse_CV,
            smooth_GCV
          )
      }
      
      # else if (method == "power" & alpha_orth == "TRUE") {
      else if (method == "eigen" ||
               alpha_orth == "TRUE") {
        # Adjust the vector to match the required lengths if they are incorrect
        if (is.vector(smooth_tuning) &
            !is.list(smooth_tuning)) {
          if (smooth_GCV == FALSE) {
            if (length(smooth_tuning) != mvmfd_obj$nvar) {
              warning(
                "The length of 'smooth_tuning' did not match number of variables and has been adjusted accordingly.",
                call. = FALSE
              )
              smooth_tuning <-
                rep(smooth_tuning, length.out = mvmfd_obj$nvar)
            }
            smooth_tuning <-
              lapply(1:mvmfd_obj$nvar, function(i)
                smooth_tuning[i])
          }
          else{
            warning(
              "The length of 'smooth_tuning' did not match number of variables and has been adjusted accordingly.",
              call. = FALSE
            )
            smooth_tuning <-
              replicate(mvmfd_obj$nvar, smooth_tuning, simplify = FALSE)
          }
        }
        
        # Adjust the matrix to match the required if they are incorrect
        else if (is.matrix(smooth_tuning)) {
          if (smooth_GCV == FALSE) {
            if (!all(dim(smooth_tuning) == c(mvmfd_obj$nvar, 1))) {
              smooth_tuning <-
                smooth_tuning[rep(1:nrow(smooth_tuning), length.out = mvmfd_obj$nvar), rep(1:ncol(smooth_tuning), length.out = 1)]
              smooth_tuning <-
                as.list(smooth_tuning)
              warning(
                "The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.",
                call. = FALSE
              )
            } else{
              smooth_tuning <- as.list(smooth_tuning)
            }
          }
          else{
            if (dim(smooth_tuning)[1] != mvmfd_obj$nvar) {
              smooth_tuning <-
                smooth_tuning[rep(1:nrow(smooth_tuning), length.out = mvmfd_obj$nvar), , drop = FALSE][1:mvmfd_obj$nvar, , drop = FALSE]
              smooth_tuning <-
                split(smooth_tuning, row(smooth_tuning))
              warning(
                "The dimensions of 'smooth_tuning' did not match the expected size and have been adjusted accordingly.",
                call. = FALSE
              )
            }
            else{
              smooth_tuning <- split(smooth_tuning, row(smooth_tuning))
            }
          }
        }
        
        # Adjust the list length and element sizes to match the required dimensions if they are incorrect
        else if (is.list(smooth_tuning)) {
          if (smooth_GCV == FALSE) {
            if (length(smooth_tuning) != mvmfd_obj$nvar) {
              warning("Adjusting 'smooth_tuning' to match 'mvmfd_obj$nvar'.",
                      call. = FALSE)
              smooth_tuning <-
                rep(smooth_tuning, length.out = mvmfd_obj$nvar)
            }
            smooth_tuning <-
              lapply(smooth_tuning, function(vec) {
                if (length(vec) != 1) {
                  warning("Adjusting vector length in 'smooth_tuning' to match 'ncomp'.",
                          call. = FALSE)
                  vec <-
                    rep(vec, length.out = 1)
                }
                vec
              })
          }
          else{
            if (length(smooth_tuning) != mvmfd_obj$nvar) {
              warning("Adjusting 'smooth_tuning' to match 'mvmfd_obj$nvar'.",
                      call. = FALSE)
              smooth_tuning <-
                rep(smooth_tuning, length.out = mvmfd_obj$nvar)[1:mvmfd_obj$nvar]
            }
          }
        }
        
        if (!is.null(smooth_tuning)) {
          names(smooth_tuning) <- paste0("var", 1:mvmfd_obj$nvar)
        }
        if (method == "power") {
          result <-
            joint_power(
              mvmfd_obj = mvmfd_obj,
              n = ncomp,
              smooth_tuning = smooth_tuning,
              centerfns = centerfns,
              alpha_orth = alpha_orth,
              smooth_tuning_type = smoothing_type
            )
        } else{
          result <-
            eigen_approach(
              mvmfd_obj = mvmfd_obj,
              n = ncomp,
              alpha = smooth_tuning,
              centerfns = centerfns,
              penalty_type = smoothing_type
            )
        }
        
        # result <- joint_power(mvmfd_obj = mvmfd_obj, n = ncomp, smooth_tuning = smooth_tuning, centerfns = centerfns, alpha_orth = alpha_orth, smooth_tuning_type = smoothing_type)
      }
      # else if (method == "eigen") {
      #   result <- eigen_approach(mvmfd_obj = mvmfd_obj, n = ncomp, alpha = smooth_tuning, centerfns = centerfns, penalty_type = smoothing_type)
      # }
      coef <- result[[1]]
      pcmfd <- list()
      for (i in 1:mvmfd_obj$nvar) {
        if (mvmfd_obj$basis$dimSupp[i] == 1) {
          pcmfd[[i]] <-
            Mfd(X = coef[[i]],
                mdbs = mvmfd_obj$basis$basis[[i]],
                method = "coefs")
        } else {
          coef_new <-
            array(coef[[i]], dim = c(mvmfd_obj$basis$basis[[i]]$nbasis, ncol(coef[[i]])))
          pcmfd[[i]] <-
            Mfd(X = coef_new,
                mdbs = mvmfd_obj$basis$basis[[i]],
                method = "coefs")
        }
      }
      out <- Mvmfd(pcmfd)
      if (mvmfd_obj$nvar == 1) {
        private$.pc_mfd <- pcmfd[[1]]
      } else {
        private$.pc_mfd <- Mvmfd(pcmfd)
      }
      private$.lsv <- result[[2]]
      private$.values <- result[[3]]
      private$.smooth_tuning <-
        result[[4]]
      if (alpha_orth == "FALSE" &&
          method == "power") {
        private$.sparse_tuning <- result[[5]]
        private$.CVs <- result[[6]]
        private$.GCVs <- result[[7]]
      } else{
        private$.GCVs <- result[[5]]
      }
      # private$.sparse_tuning <- result[[5]]
      private$.mean_mfd <-
        mean(mvmfd_obj)
    }
  ),
  active = list(
    pc_mfd = function(value) {
      if (missing(value)) {
        private$.pc_mfd
      } else {
        stop("`$pc_mfd` is read only", call. = FALSE)
      }
    },
    lsv = function(value) {
      if (missing(value)) {
        private$.lsv
      } else {
        stop("`$lsv` is read only", call. = FALSE)
      }
    },
    values = function(value) {
      if (missing(value)) {
        private$.values
      } else {
        stop("`$sigma` is read only", call. = FALSE)
      }
    },
    smooth_tuning = function(value) {
      if (missing(value)) {
        private$.smooth_tuning
      } else {
        stop("`$smooth_tuning` is read only", call. = FALSE)
      }
    },
    sparse_tuning = function(value) {
      if (missing(value)) {
        private$.sparse_tuning
      } else {
        stop("`$sparse_tuning` is read only", call. = FALSE)
      }
    },
    GCVs = function(value) {
      if (missing(value)) {
        private$.GCVs
      } else {
        stop("`$GCVs` is read only", call. = FALSE)
      }
    },
    CVs = function(value) {
      if (missing(value)) {
        private$.CVs
      } else {
        stop("`$CVs` is read only", call. = FALSE)
      }
    },
    mean_mfd = function(value) {
      if (missing(value)) {
        private$.mean_mfd
      } else {
        stop("`$mean_mfd` is read only", call. = FALSE)
      }
    }
  ),
  private = list(
    .pc_mfd = NULL,
    .lsv = NULL,
    .values = NULL,
    .smooth_tuning = NULL,
    .sparse_tuning = NULL,
    .GCVs = NULL,
    .CVs = NULL,
    .mean_mfd = NULL
  )
)

#' @rdname remfpca
#' @seealso \code{\link{mvmfd}}

#' @title A Class for 'remfpca' objects
#'
#' @description
#' The `remfpca` class represents regularized functional principal components ('ReMFPCs') components.
#'
#' @param mvmfd_obj An `mvmfd` object representing the multivariate functional data.
#' @param method A character string specifying the approach to be used for MFPCA computation.
#'               Options are "power" (the default), which uses the power algorithm, or "eigen",
#'               which uses the eigen decomposition approach.
#' @param ncomp The number of functional principal components to retain.
#' @param smooth_tuning A list or vector specifying the smoothing regularization parameter(s) for each variable.
#'                      If NULL, non-smoothing MFPCA is estimated.
#' @param sparse_tuning A list or vector specifying the sparsity regularization parameter(s) for each variable.
#'                      If NULL, non-sparse MFPCA is estimated.
#' @param centerfns Logical indicating whether to center the functional data before analysis. Default is TRUE.
#' @param alpha_orth Logical indicating whether to perform orthogonalization of the regularization parameters.
#'                   If `method` is "power", setting `alpha_orth = FALSE` (default) uses the sequential power approach,
#'                   while setting `alpha_orth = TRUE` uses the joint power approach.
#' @param smoothing_type The type of smoothing penalty to be applied on the estimated functional PCs. The types "basispen" and "coefpen" is supported. Default is "basispen".
#' @param sparse_type The type of sparse penalty to be applied on the estimated functional PCs. The types "soft-threshold", "hard-threshold" and "SCAD" is supported. Default is "soft-threshold".
#' @param K_fold  An integer specifying the number of folds in the sparse cross-validation process. Default is 30.
#' @param sparse_CV @param sparse_CV Logical indicating whether cross-validation should be applied to select the optimal sparse tuning parameter in sequential power approach.
#'                                        If `sparse_CV = TRUE`, a series of tuning parameters should be provided as a vector with positive number with max equals to number of subjects.
#'                                        If `sparse_CV = FALSE`, specific tuning parameters are given directly to each principal components. Tuning parameters should be provided as a vector with length equal to `ncomp`.
#'                                        If the dimensions of input tuning parameters are incorrect, it will be converted to a list internally, and a warning will be issued.
#' @param smooth_GCV @param smooth_GCV Logical indicating whether generalized cross-validation should be applied to select the optimal smooth tuning parameter.
#'                                        If `smooth_GCV = TRUE`, a series of tuning parameters should be provided as a list with length equal to the number of variables.
#'                                        If a list with incorrect dimensions is provided, it will be converted to a correct list internally, and a warning will be issued.
#'                                        If `smooth_GCV = FALSE`, specific tuning parameters are given directly. If `method` is "power" and `alpha_orth = FALSE` (sequential power),
#'                                        tuning parameters should be provided as a list with length equal to the number of variables, where each element is a vector of length `ncomp`.
#'                                        If `method` is "power" and `alpha_orth = TRUE` (joint power), tuning parameters should be provided as a vector with length equal to the number of variables.
#'                                        If the dimensions of input tuning parameters are incorrect, it will be converted to a list internally, and a warning will be issued.
#' @export
Remfpca <-
  function(mvmfd_obj,
           method = "power",
           ncomp,
           smooth_tuning = NULL,
           sparse_tuning = NULL,
           centerfns = TRUE,
           alpha_orth = FALSE,
           smoothing_type = "basispen",
           sparse_type = "soft",
           K_fold = 30,
           sparse_CV = TRUE,
           smooth_GCV = TRUE) {
    remfpca$new(
      mvmfd_obj,
      method,
      ncomp,
      smooth_tuning,
      sparse_tuning,
      centerfns,
      alpha_orth,
      smoothing_type,
      sparse_type,
      K_fold,
      sparse_CV,
      smooth_GCV
    )
  }
#' @rdname remfpca
