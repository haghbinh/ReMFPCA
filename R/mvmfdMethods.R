#' Addition of two mvmfd objects
#'
#' @export
#'
#' @param obj1 An mvmfd object
#' @param obj2 An optional mvmfd object
#' @return An mvmfd object
#'
`+.mvmfd` <- function(obj1, obj2 = NULL) {
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

#' Subtraction of two mvmfd objects
#'
#' @export
#'
#' @param obj1 An mvmfd object
#' @param obj2 An optional mvmfd object
#' @return An mvmfd object
#'
`-.mvmfd` <- function(obj1, obj2 = NULL) {
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

#' Multiplication of an mvmfd object with a scalar
#'
#' @export
#'
#' @param obj1 An mvmfd object or a scalar
#' @param obj2 An mvmfd object or a scalar
#' @return An mvmfd object
#'
`*.mvmfd` <- function(obj1, obj2) {
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

#' Plotting method for mvmfd objects
#'
#' @export
#'
#' @param mvmfd_obj An mvmfd object
#' @param xlab Label for the x-axis
#' @param ylab Label for the y-axis
#' @param ... Additional arguments for the plot function
#'
plot.mvmfd <- function(mvmfd_obj, xlab = NULL, ylab = NULL, ...) {
  old <- par(no.readonly = TRUE, mfrow = c(1, 1))
  p <- mvmfd_obj$nvar
  par(mfrow = c(p, 1))
  if (is.null(ylab)) ylab <- paste("Variable ", 1:p)
  if (is.null(xlab)) xlab <- rep("time", p)
  for (i in 1:p) {
    plot(mvmfd_obj[, i], ylab = ylab[i], xlab = xlab[i], ...)
  }
  par(mfrow = c(1, 1))
  on.exit(options(old))
}

#' Bivariate plot for mvmfd objects
#'
#' @export
#'
#' @param mvmfd_obj An mvmfd object
#' @param type Type of plot ('l' for lines, 'p' for points, etc.)
#' @param lty Line type
#' @param xlab Label for the x-axis
#' @param ylab Label for the y-axis
#' @param main Main title
#' @param ... Additional arguments for the matplot function
#'
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

#' Mean of each variable in an mvmfd object
#'
#' @export
#'
#' @param mvmfd_obj An mvmfd object
#' @return An mvmfd object
#'
mean.mvmfd <- function(mvmfd_obj) {
  p <- mvmfd_obj$nvar
  mvlist <- lapply(1:p, function(j) mean(mvmfd_obj[, j]))
  return(Mvmfd(mvlist))
}

#' Inner product of two mvmfd objects
#'
#' @export
#'
#' @param mvmfd_obj1 An mvmfd object
#' @param mvmfd_obj2 An mvmfd object
#' @return A matrix of inner products
#'
inprod_mvmfd <- function(mvmfd_obj1, mvmfd_obj2) {
  p <- mvmfd_obj1$nvar
  if (p != mvmfd_obj2$nvar) stop("The number of variables must be equal.")
  m <- mvmfd_obj1$nobs
  n <- mvmfd_obj2$nobs
  inpr <- matrix(0, nrow = m, ncol = n)
  for (j in 1:p) {
    inpr <- inpr + inprod_mfd(mvmfd_obj1[, j], mvmfd_obj2[, j])
  }
  return(inpr)
}

#' Norm of an mvmfd object
#'
#' @export
#'
#' @param mvmfd_obj An mvmfd object
#' @return The norm of the mvmfd object
#'
norm_mvmfd <- function(mvmfd_obj) {
  return(as.numeric(sqrt(inprod_mvmfd(mvmfd_obj, mvmfd_obj))))
}
