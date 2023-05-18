#' A Class for ReMFPCA objects
#' @description The \code{remfpca} class represents functional data ...
#' @field basis a basismfd object
#' @field coeff a matrix with nrow=subjects and ncol=total number of basis ...
#'
#' @examples
#' x <- 1
#'
#' @importFrom fda getbasispenalty
#' @importFrom Matrix bdiag
#' @importFrom expm sqrtm
#' @export
remfpca <- R6::R6Class("remfpca",
  public = list(
    initialize = function(mvmfd_obj, method = "eigen", ncomp, alpha = NULL, centerfns = TRUE, alpha_orth = TRUE, lambda_type = "variable",penalty_type = "coefpen") {
      if (is.numeric(alpha)) alpha <- as.list(alpha)
      if (is.mfd(mvmfd_obj)) mvmfd_obj <- mvmfd$new(mvmfd_obj)
      if (method == "power") {
        result <- power_algo_fun(mvmfd_obj = mvmfd_obj, n = ncomp, alpha = alpha, centerfns = centerfns, alpha_orth = alpha_orth, lambda_type = lambda_type,penalty_type = penalty_type)
      } else if (method == "halfsmooth") {
        result <- half_smoothing_approach(mvmfd_obj = mvmfd_obj, n = ncomp, alpha = alpha, centerfns = centerfns)
      } else if (method == "eigen") {
        result <- eigen_approach(mvmfd_obj = mvmfd_obj, n = ncomp, alpha = alpha, centerfns = centerfns,penalty_type = penalty_type)
      }
      coef <- result[[1]]
      pcmfd <- list()
      for (i in 1:mvmfd_obj$nvar) {
        if (mvmfd_obj$basis$dimSupp[i] == 1) {
          pcmfd[[i]] <- Mfd(X = coef[[i]], mdbs = mvmfd_obj$basis$basis[[i]], method = "coefs")
        } else {
          coef_new <- array(coef[[i]], dim = c(mvmfd_obj$basis$basis[[i]]$nbasis, ncol(coef[[i]])))
          pcmfd[[i]] <- Mfd(X = coef_new, mdbs = mvmfd_obj$basis$basis[[i]], method = "coefs")
        }
      }
      out <- Mvmfd(pcmfd)
      if(mvmfd_obj$nvar==1) {
        private$.pc_mfd <- pcmfd[[1]]
        } else {
        private$.pc_mfd <- Mvmfd(pcmfd) 
      }
      private$.lsv <- result[[2]]
      private$.values <- result[[3]]
      private$.alpha <- result[[4]]
      private$.GCVs <- result[[5]]
      private$.mean_mfd <- mean(mvmfd_obj)
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
    alpha = function(value) {
      if (missing(value)) {
        private$.alpha
      } else {
        stop("`$alpha` is read only", call. = FALSE)
      }
    },
    GCVs = function(value) {
      if (missing(value)) {
        private$.GCVs
      } else {
        stop("`$GCVs` is read only", call. = FALSE)
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
    .alpha = NULL,
    .GCVs = NULL,
    .mean_mfd = NULL
  )
)

#' @export
Remfpca <- function(mvmfd_obj, method = "eigen", ncomp, alpha = NULL, centerfns = TRUE, alpha_orth = TRUE, lambda_type = "variable",penalty_type = "coefpen") remfpca$new(mvmfd_obj, method, ncomp, alpha, centerfns, alpha_orth, lambda_type,penalty_type)
