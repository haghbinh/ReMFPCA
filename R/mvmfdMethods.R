
#' @export
"+.mvmfd" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
    return(obj1)
  }
  p <- obj1$nvar
  mvlist <- list()
  for (j in 1:p) {
    mvlist[[j]] <- obj1[, j] + obj2[, j]
  }
  return(Mvmfd(mvlist))
}

#' @export
"-.mvmfd" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
    return((-1) * obj1)
  }
  p <- obj1$nvar
  mvlist <- list()
  for (j in 1:p) {
    mvlist[[j]] <- obj1[, j] + (-1) * obj2[, j]
  }
  return(Mvmfd(mvlist))
}

#' @export
"*.mvmfd" <- function(obj1, obj2) {
  if (xor(is.mvmfd(obj1), is.mvmfd(obj2))) {
    if (xor(is.double(obj1), is.double(obj2))) {
      if (is.double(obj1)) {
        temp <- obj1
        obj1 <- obj2
        obj2 <- temp
      }
    }
    p <- obj1$nvar
    mvlist <- list()
    for (j in 1:p) {
      mvlist[[j]] <- obj1[, j] * obj2
    }
  } else {
    stop("One object must be an mvmfd, and the other one a scalar")
  }
  return(Mvmfd(mvlist))
}


#' @export
plot.mvmfd <- function(mvfd_obj, type = "", xlab = "", main = "", ...) {
  p <- mvmfd_obj$nvar
  par(mfrow = c(p, 1))
  for (i in 1:p) {
    plot(mvfd_obj[, i], ylab = paste("Variable ", i), obs, xlab, main, ...)
  }
}

#' @export
bimfdplot <- function(mvmfd_obj, type = "l", lty = 1, xlab = "", ylab = "", main = "", ...) {
  nvar <- mvmfd_obj$nvar
  stopifnot(nvar == 2)
  stopifnot(all(mvmfd_obj$basis$supp[[1]] == mvmfd_obj$basis$supp[[2]]))
  supp <- mvmfd_obj$basis$supp[[1]]
  x_grids <- seq(supp[1, 1], supp[2, 1], len = 1000)
  X <- mvmfd_obj[, 1]$eval(x_grids)
  Y <- mvmfd_obj[, 2]$eval(x_grids)
  matplot(X, Y, type = type, lty = lty, xlab = xlab, ylab = ylab, main = main, ...)
}


#' @export
mean.mvmfd <- function(mvmfd_obj) {
  p <- mvmfd_obj$nvar
  mvlist <- lapply(1:p, function(j) mean(mvmfd_obj[, j]))
  return(Mvmfd(mvlist))
}

#' @export
inprod_mvmfd <- function(mvmfd_obj1, mvmfd_obj2) {
  p <- mvmfd_obj1$nvar
  if (p != mvmfd_obj2$nvar) stop("The number of variablles must be equal.")
  m <- mvmfd_obj1$nobs
  n <- mvmfd_obj2$nobs
  inpr <- matrix(0, nrow = m, ncol = n)
  for (j in 1:p) {
    inpr <- inpr + inprod_mfd(mvmfd_obj1[, j], mvmfd_obj2[, j])
  }
  return(inpr)
}

#' @export
norm_mvmfd <- function(mvmfd_obj) {
  return(as.numeric(sqrt(inprod_mvmfd(mvmfd_obj, mvmfd_obj))))
}
