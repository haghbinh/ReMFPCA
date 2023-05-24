#' @title Length of an 'mfd' object
#'
#' @param mfd_obj An 'mfd' object
#' @return The number of observations in the 'mfd' object
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
length.mfd <- function(mfd_obj) {
  return(mfd_obj$nobs)
}

#' @title Plot an 'mfd' object
#'
#' @param mfd_obj An 'mfd' object
#' @param obs Observation number to plot (default: 1)
#' @param xlab x-axis label (default: "")
#' @param ylab y-axis label (default: "")
#' @param main Title of the plot (default: "")
#' @param type Type of plot (default: "l")
#' @param lty Line type (default: 1)
#' @param ... Additional arguments passed to the plot function
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
plot.mfd <- function(mfd_obj, obs = 1, xlab = "", ylab = "", main = "", type = "l", lty = 1, ...) {
  dimSupp <- mfd_obj$basis$dimSupp
  supp <- mfd_obj$basis$supp
  x_grids <- seq(supp[1, 1], supp[2, 1], len = 1000)
  if (dimSupp == 1) {
    X <- mfd_obj$eval(x_grids)
    matplot(x_grids, X, type = type, lty = lty, xlab = xlab, ylab = ylab, main = main, ...)
  } else {
    y_grids <- seq(supp[1, 2], supp[2, 2], len = 100)
    X <- mfd_obj$eval(list(x_grids, y_grids))[, , obs]
    image(X, xlab = xlab, ylab = ylab, axes = FALSE, main = paste(main, " Observation:", obs))
    axis(side = 1, at = seq(0, 1, len = 10), labels = round(seq(supp[1, 1], supp[2, 1], len = 10), 1))
    axis(side = 2, at = seq(0, 1, len = 10), labels = round(seq(supp[1, 2], supp[2, 2], len = 10), 1))
  }
}

#' @title Compute the mean of an 'mfd' object
#'
#' @param mfd_obj An 'mfd' object
#' @return An 'mfd' object representing the mean
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
mean.mfd <- function(mfd_obj) {
  cof <- apply(mfd_obj$coefs, 1, mean)
  bs <- mfd_obj$basis
  return(Mfd(X = cof, mdbs = bs, method = "coefs"))
}

#' @title Compute the standard deviation of an 'mfd' object
#'
#' @param mfd_obj An 'mfd' object
#' @return An 'mfd' object representing the standard deviation
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @importFrom stats sd
#' @export
sd.mfd <- function(mfd_obj) {
  cof <- apply(mfd_obj$coefs, 1, sd)
  bs <- mfd_obj$basis
  return(Mfd(X = cof, mdbs = bs, method = "coefs"))
}

#' @title Compute the inner product between two 'mfd' objects
#'
#' @param mfd_obj1 An 'mfd' object
#' @param mfd_obj2 An 'mfd' object
#' @return The inner product between the two 'mfd' objects
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @importFrom fda fd inprod
#' @export
inprod_mfd <- function(mfd_obj1, mfd_obj2) {
  bs1 <- mfd_obj1$basis$basis[[1]]
  bs2 <- mfd_obj2$basis$basis[[1]]
  cof1 <- mfd_obj1$coefs
  cof2 <- mfd_obj2$coefs
  fd1 <- fd(coef = cof1, basisobj = bs1)
  fd2 <- fd(coef = cof2, basisobj = bs2)
  inpr <- fda::inprod(fd1, fd2)
  return(inpr)
}

#' @title Compute the norm of an 'mfd' object
#'
#' @param mfd_obj An 'mfd' object
#' @return The norm of the 'mfd' object
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
norm_mfd <- function(mfd_obj) {
  return(as.numeric(sqrt(inprod_mfd(mfd_obj, mfd_obj))))
}

#' @title Add two 'mfd' objects
#'
#' @param obj1 An 'mfd' object
#' @param obj2 An 'mfd' object or a scalar
#' @return The sum of the two 'mfd' objects
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
"+.mfd" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
    return(obj1)
  }
  if (xor(is.mfd(obj1), is.mfd(obj2))) {
    if (!xor(is.double(obj1), is.double(obj2))) stop("At least one object must be an mfd, and the other one can be a scalar")
    if (is.double(obj1)) {
      temp <- obj1
      obj1 <- obj2
      obj2 <- temp
    }
    coef <- obj1$coefs + obj2
  } else {
    if (length(obj1) == length(obj2)) {
      if (obj1$basis$dimSupp != obj2$basis$dimSupp) stop("Two objects must have the same basis dimSupp.")
      if (obj1$basis$dimSupp == 1) {
        if (obj1$basis$basis[[1]]$type != obj2$basis$basis[[1]]$type) stop("Two objects must have the same basis types.")
      } else {
        for (i in 1:obj1$basis$dimSupp) {
          if (obj1$basis$basis[[i]]$type != obj2$basis$basis[[i]]$type) stop("Two objects must have the same basis types.")
        }
      }
      coef <- obj1$coefs + obj2$coefs
    } else {
      if (length(obj1) > 1 & length(obj2) > 1) stop("Two objects must have the same basis dimSupp.")
      if (length(obj1) == 1 & length(obj2) > 1) {
        temp <- obj1
        obj1 <- obj2
        obj2 <- temp
      }
      nobs <- obj1$nobs
      I_n <- matrix(1L, ncol = nobs, nrow = 1)
      coef <- obj1$coefs + as.matrix(obj2$coefs) %*% I_n
    }
  }
  return(mfd$new(X = coef, mdbs = obj1$basis$clone(), method = "coefs"))
}

#' @title Multiply two 'mfd' objects
#'
#' @param obj1 An 'mfd' object
#' @param obj2 An 'mfd' object or a scalar
#' @return The product of the two 'mfd' objects
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
"*.mfd" <- function(obj1, obj2) {
  if (xor(is.mfd(obj1), is.mfd(obj2))) {
    if (xor(is.double(obj1), is.double(obj2))) {
      if (is.double(obj1)) {
        temp <- obj1
        obj1 <- obj2
        obj2 <- temp
      }
    }
    coef <- obj1$coefs * obj2
  } else {
    stop("One object must be an mfd, and the other one a scalar")
  }
  return(mfd$new(X = coef, mdbs = obj1$basis$clone(), method = "coefs"))
}

#' @title Subtract two 'mfd' objects
#'
#' @param obj1 An 'mfd' object
#' @param obj2 An 'mfd' object or a scalar
#' @return The difference between the two 'mfd' objects
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
"-.mfd" <- function(obj1, obj2 = NULL) {
  if (is.null(obj2)) {
    return((-1) * obj1)
  }
  if (xor(is.mfd(obj1), is.mfd(obj2))) {
    if (!xor(is.double(obj1), is.double(obj2))) stop("At least one object must be an mfd, and the other one can be a scalar")
    if (is.double(obj1)) {
      coef <- obj1 + (-1) * obj2$coefs
      obj1 <- obj2
    } else {
      coef <- obj1$coefs + (-obj2)
    }
  } else {
    return(obj1 + (-1) * obj2)
  }
  return(mfd$new(X = coef, mdbs = obj1$basis$clone(), method = "coefs"))
}

#' @title Extract subsets of an 'mfd' object
#'
#' @param mfd_obj An 'mfd' object
#' @param i An index or indices specifying the subsets to extract
#' @return An 'mfd' object containing the specified subsets
#' @seealso \code{\link{basismfd}}, \code{\link{mfd}}
#' @export
"[.mfd" <- function(mfd_obj, i = "index") {
  if(is.null(i)) i <- 1:mfd_obj$nobs
  if (max(i) > mfd_obj$nobs | min(i) < 1) stop(" subscript i out of bounds")
  bs <- mfd_obj$basis$clone()
  if(mfd_obj$basis$dimSupp==1){
    coef <- mfd_obj$coefs[, i]
  }else{
    coef <- mfd_obj$coefs[, , i]
  }
  return(mfd$new(X = coef, mdbs = bs, method = "coefs"))
}
