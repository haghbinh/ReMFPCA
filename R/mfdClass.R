#' Multivariate Functional Data object Class
#'
#' This function is used to create a multivariate functional data (mfd) objects from lists of discretely sampled data,
#'  basis specifications, and grid elements which provide the domain that each variable is observed over. Each variable is assumed to be observed over a regular and equidistant grid.
#'  In addition, each variable in the mfd is assumed to be observed over a one or two-dimensional domain.
#' @param X  A list of length \code{p} where \code{p} is a positive integer specifying the number of variables observed in the mfd. Each entry in the list should be a matrix or
#' an array. For a variable observed over a one-dimensional domain, the list entry should be an \code{m} by \code{N} matrix where \code{m} is the number of sampling points for each observation and \code{N} is the number of observations.
#' For a variable observed over a two-dimensional domain, the list entry should be an \code{m_1} by \code{m_2} by \code{N} array where \code{m_1} is the number of sampling points in the horizontal direction and
#' \code{m_2} is the number of sampling points in the vertical direction.
#' @param B A list of length \code{p}. Each entry in the list should be either a matrix specifying the basis for each variable or each list entry should be a list specifying the number of basis elements and desired basis type to be used in the smoothing process.
#' For a variable observed over a one-dimensional domain, the list entry should be a \code{m} by \code{d} matrix or the list entry should be a list of the form, \code{list(d, "basistype")}, where \code{d} specifies the number basis elements and \code{"basistype"} is a string that is set to \code{"bspline"} or \code{"fourier"} to specify the type of basis.
#' For a variable observed over a two-dimensional domain, the list entry should be a \code{(m_1)(m_2)} by \code{(d_1)(d_2)} matrix or the list entry should be a list of the form,
#' \code{list(d_1,d_2,"basistype1","basistype2")}, where \code{d_1} and \code{d_2} specify the number of basis elements in the horizontal and vertical directions respectively, and \code{"basistype1"} and \code{"basistype2"} are set to be \code{"bspline"} or \code{"fourier"} which specifies the type of basis in the horizontal direction, and \code{"basistype2"} is set to be \code{"bspline"} or \code{"fourier"} to specify the type of basis in the horizontal and vertical directions respectively.
#' @param grid A list of length \code{p}. Each entry in the list should either be a numeric or a list of numeric elements depending on the dimension of the domain the variable is observed over.
#' For a variable observed over a one-dimensional domain, the list entry should be either an ordered (from smallest to largest) numeric of length \code{m} giving the sampling points or a numeric of the form, \code{c(x_1,x_2)} where \code{x_1} and \code{x_2} are the smallest and largest values attained in the domain of the variable respectively. In addition, these list entries can also be provided in a list.
#' For a variable observed over a two-dimensional domain, the list entry should be either a list of the form, \code{list(u,v)}, or of the form, \code{list(c(x_1,x_2),v)}, or of the form, \code{list(u,c(x_3,x_4))}, or of the form, \code{list(c(x_1,x_2),c(x_3,x_4))} where \code{u} is a numeric with minimium and maximum values of \code{x_1} and \code{x_2} respectively and \code{v} is a numeric with minimium and maximum values of \code{x_3} and \code{x_4} respectively.
#' @importFrom fda fd inprod eval.fd smooth.basis is.fd create.bspline.basis create.fourier.basis eval.basis
#' @importFrom methods is new
#' @example 
#' 
#' @export
mfd <- function(X, B, grid) {
  mfd <- new(Class = "mfd", C = X, B = B, grid = grid)
  p <- length(X)
  # Loop over each variable in the mfd
  for (j in 1L:p) {
    # Deciding whether or not the variable corresponds to an mfd observed over a two-dimensional domain or a one-dimensional domain.
    if (class(grid[[j]]) == "list" && length(grid[[j]]) == 2) {
      if (length(grid[[j]][[1]]) == 2 && length(grid[[j]][[2]]) == 2) {
        # If the grid is specified with min and max values in both horizontal and vertical axes.
        # Example: list(c(1,2),c(3,4))
        
        u <- seq(grid[[j]][[1]][1], grid[[j]][[1]][2], length.out = ncol(X[[j]]))
        v <- seq(grid[[j]][[2]][1], grid[[j]][[2]][2], length.out = nrow(X[[j]]))
      } else if (length(grid[[j]][[1]]) == 2 && length(grid[[j]][[2]]) > 2) {
        # If the grid is specified with min and max values in the horizontal axis and with the actual grid in the vertical axis.
        # Example: list(c(1,2),v)
        
        u <- seq(grid[[j]][[1]][1], grid[[j]][[1]][2], length.out = ncol(X[[j]]))
        v <- grid[[j]][[2]]
      } else if (length(grid[[j]][[1]]) > 2 && length(grid[[j]][[2]]) == 2) {
        # This case is similar as seen in the previous else if statement
        # Example: list(u,c(1,2))
        
        u <- grid[[j]][[1]]
        v <- seq(grid[[j]][[2]][1], grid[[j]][[2]][2], length.out = nrow(X[[j]]))
      } else {
        
        # If the grid is specified with the actual grid in the horizontal and vertical axes.
        # Example: list(u,v)
        u <- seq(min(grid[[j]][[1]]), max(grid[[j]][[1]]), length.out = length(grid[[j]][[1]]))
        v <- seq(min(grid[[j]][[2]]), max(grid[[j]][[2]]), length.out = length(grid[[j]][[2]]))
      }
      M_x <- length(u)
      M_y <- length(v)
      
      # Create the grid using u and v
      mfd@grid[[j]] <- cbind(rep(u, each = length(v)), rep(v, length(u)))
      
      # Reshape mfd of Images from matricies to vectors
      if (is.matrix(X[[j]])) X[[j]] <- array(X[[j]], dim = c(M_x, M_y, 1))
      X[[j]] <- matrix(aperm(X[[j]], c(2, 1, 3)), nrow = M_x * M_y)
      
      # If the variable is observed over a one-dimensional domain
    } else if (class(grid[[j]]) == "list" && length(grid[[j]]) == 1 && length(grid[[j]][[1]]) > 2) {
      # If the grid is supplied in a list, the variable corresponds with a one-dimensional domain, and the list element is the grid.
      # Example: list(u)
      M_x <- length(grid[[j]][[1]])
      u <- seq(min(grid[[j]][[1]]), max(grid[[j]][[1]]), length.out = M_x)
      mfd@grid[[j]] <- matrix(u, ncol = 1)
      X[[j]] <- as.matrix(X[[j]])
    } else if (class(grid[[j]]) == "list" && length(grid[[j]]) == 1 && length(grid[[j]][[1]]) == 2) {
      
      # If the grid is supplied in a list, the variable corresponds with a one-dimensional domain, and the list element is the min and max values of the grid.
      # Example: list(c(1,2))
      u <- seq(grid[[j]][[1]][1], grid[[j]][[1]][2], length.out = nrow(X[[j]]))
      mfd@grid[[j]] <- matrix(u, ncol = 1)
      X[[j]] <- as.matrix(X[[j]])
    } else if (class(grid[[j]]) == "numeric" && length(grid[[j]]) == 2) {
      
      # If the grid specified with the min and max values of the grid.
      # Example: c(1,2)
      u <- seq(grid[[j]][1], grid[[j]][2], length.out = nrow(X[[j]]))
      mfd@grid[[j]] <- matrix(u, ncol = 1)
      X[[j]] <- as.matrix(X[[j]])
    } else {
      # If the grid specified with the grid itself.
      # Example: u
      mfd@grid[[j]] <- matrix(grid[[j]], ncol = 1)
      X[[j]] <- as.matrix(X[[j]])
    }
    
    
    
    # Generating basis matrices
    if (class(B[[j]]) == "list" && B[[j]][[2]] == "bspline" && ncol(mfd@grid[[j]]) == 1) {
      # Generating a bspline basis for variables whose domain is one-dimensional using a supplied list.
      # Example: list(10,"bspline")
      mfd@B[[j]] <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]]), max(mfd@grid[[j]]), length.out = nrow(mfd@grid[[j]])), basisobj = fda::create.bspline.basis(rangeval = c(min(mfd@grid[[j]]), max(mfd@grid[[j]])), nbasis = B[[j]][[1]]))
      B[[j]] <- mfd@B[[j]]
    } else if (class(B[[j]]) == "list" && B[[j]][[2]] == "fourier" && ncol(mfd@grid[[j]]) == 1) {
      # Generating a fourier basis for variables whose domain is one-dimensional using a supplied list.
      # Example: list(10,"fourier")
      mfd@B[[j]] <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]]), max(mfd@grid[[j]]), length.out = nrow(mfd@grid[[j]])), basisobj = fda::create.fourier.basis(rangeval = c(min(mfd@grid[[j]]), max(mfd@grid[[j]])), nbasis = B[[j]][[1]]))
      B[[j]] <- mfd@B[[j]]
    }
    
    if (class(B[[j]]) == "list" && B[[j]][[3]] == "bspline" && B[[j]][[4]] == "bspline" && length(grid[[j]]) == 2) {
      # Generating a bspline basis for variables whose domain is two-dimensional using a supplied list.
      # Example: list(list(10,"bspline"),list(12,"bspline"))
      b_1 <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]][, 1]), max(mfd@grid[[j]][, 1]), length.out = M_x), basisobj = fda::create.bspline.basis(rangeval = c(min(mfd@grid[[j]][, 1]), max(mfd@grid[[j]][, 1])), nbasis = B[[j]][[1]]))
      b_2 <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]][, 2]), max(mfd@grid[[j]][, 2]), length.out = M_y), basisobj = fda::create.bspline.basis(rangeval = c(min(mfd@grid[[j]][, 2]), max(mfd@grid[[j]][, 2])), nbasis = B[[j]][[2]]))
      mfd@B[[j]] <- kronecker(b_1, b_2)
      B[[j]] <- mfd@B[[j]]
    } else if (class(B[[j]]) == "list" && B[[j]][[3]] == "bspline" && B[[j]][[4]] == "fourier" && length(grid[[j]]) == 2) {
      # Generating a bspline/fourier basis for variables whose domain is two-dimensional using a supplied list.
      # Example: list(list(10,"bspline"),list(12,"fourier"))
      b_1 <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]][, 1]), max(mfd@grid[[j]][, 1]), length.out = M_x), basisobj = fda::create.bspline.basis(rangeval = c(min(mfd@grid[[j]][, 1]), max(mfd@grid[[j]][, 1])), nbasis = B[[j]][[1]]))
      b_2 <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]][, 2]), max(mfd@grid[[j]][, 2]), length.out = M_y), basisobj = fda::create.fourier.basis(rangeval = c(min(mfd@grid[[j]][, 2]), max(mfd@grid[[j]][, 2])), nbasis = B[[j]][[2]]))
      mfd@B[[j]] <- kronecker(b_1, b_2)
      B[[j]] <- mfd@B[[j]]
    } else if (class(B[[j]]) == "list" && B[[j]][[3]] == "fourier" && B[[j]][[4]] == "bspline" && length(grid[[j]]) == 2) {
      # This case is similar to the previous else if statement
      # Example: list(list(10,"fourier"),list(12,"bspline"))
      b_1 <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]][, 1]), max(mfd@grid[[j]][, 1]), length.out = M_x), basisobj = fda::create.fourier.basis(rangeval = c(min(mfd@grid[[j]][, 1]), max(mfd@grid[[j]][, 1])), nbasis = B[[j]][[1]]))
      b_2 <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]][, 2]), max(mfd@grid[[j]][, 2]), length.out = M_y), basisobj = fda::create.bspline.basis(rangeval = c(min(mfd@grid[[j]][, 2]), max(mfd@grid[[j]][, 2])), nbasis = B[[j]][[2]]))
      mfd@B[[j]] <- kronecker(b_1, b_2)
      B[[j]] <- mfd@B[[j]]
    } else if (class(B[[j]]) == "list" && B[[j]][[3]] == "fourier" && B[[j]][[4]] == "fourier" && length(grid[[j]]) == 2) {
      # Generating a fourier basis for variables whose domain is two-dimensional using a supplied list.
      # Example: list(list(10,"fourier"),list(12,"fourier"))
      b_1 <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]][, 1]), max(mfd@grid[[j]][, 1]), length.out = M_x), basisobj = fda::create.fourier.basis(rangeval = c(min(mfd@grid[[j]][, 1]), max(mfd@grid[[j]][, 1])), nbasis = B[[j]][[1]]))
      b_2 <- fda::eval.basis(evalarg = seq(min(mfd@grid[[j]][, 2]), max(mfd@grid[[j]][, 2]), length.out = M_y), basisobj = fda::create.fourier.basis(rangeval = c(min(mfd@grid[[j]][, 2]), max(mfd@grid[[j]][, 2])), nbasis = B[[j]][[2]]))
      mfd@B[[j]] <- kronecker(b_1, b_2)
      B[[j]] <- mfd@B[[j]]
    }
    # Estimate the coefficients of each mfd variable using basis elements and observed data.
    
    mfd@C[[j]] <- solve(t(B[[j]]) %*% B[[j]]) %*% t(B[[j]]) %*% X[[j]]
    
    if (is.null(colnames(mfd@C[[j]])) == TRUE) {
      colnames(mfd@C[[j]]) <- as.character(1:ncol(mfd@C[[j]]))
    }
    
    if (is.null(rownames(mfd@grid[[j]])) == TRUE) {
      rownames(mfd@grid[[j]]) <- as.character(1:nrow(mfd@grid[[j]]))
    }
  }
  
  
  return(mfd)
}




# Check function (Jordan's code)
check_mfd <- function(object) {
  errors <- character()
  msg <- NULL
  
  # check if mfd lengths are equal
  if (length(object@C) > 1) {
    for (j_1 in 1:length(object@C)) {
      for (j_2 in j_1:length(object@C)) {
        if (length(dim(object@C[[j_1]])) == 2 && length(dim(object@C[[j_2]])) == 2) {
          if (ncol(object@C[[j_1]]) != ncol(object@C[[j_2]])) {
            msg <- paste("The length of variable", as.character(j_1), "and the length of variable", as.character(j_2), "are not equal.")
            
            errors <- c(errors, msg)
            
            msg <- NULL
          }
        } else if (length(dim(object@C[[j_1]])) == 2 && length(dim(object@C[[j_2]])) == 3) {
          if (ncol(object@C[[j_1]]) != dim(object@C[[j_2]])[3]) {
            msg <- paste("The length of variable", as.character(j_1), "and the length of variable", as.character(j_2), "are not equal.")
            
            errors <- c(errors, msg)
            
            msg <- NULL
          }
        } else if (length(dim(object@C[[j_1]])) == 3 && length(dim(object@C[[j_2]])) == 2) {
          if (dim(object@C[[j_1]])[3] != ncol(object@C[[j_2]])) {
            msg <- paste("The length of variable", as.character(j_1), "and the length of variable", as.character(j_2), "are not equal.")
            
            errors <- c(errors, msg)
            
            msg <- NULL
          }
        } else if (length(dim(object@C[[j_1]])) == 3 && length(dim(object@C[[j_2]])) == 3) {
          if (dim(object@C[[j_1]])[3] != dim(object@C[[j_2]])[3]) {
            msg <- paste("The length of variable", as.character(j_1), "and the length of variable", as.character(j_2), "are not equal.")
            
            errors <- c(errors, msg)
            
            msg <- NULL
          }
        }
      }
    }
  }
  
  # check if provided lists have same length
  if (length(object@C) != length(object@B)) {
    msg <- paste("The length of X is", as.character(length(object@C)), "and the length of B is", as.character(length(object@B)), "whereas the length of these lists should be equal.")
    
    errors <- c(errors, msg)
    
    msg <- NULL
  }
  if (length(object@C) != length(object@grid)) {
    msg <- paste("The length of X is", as.character(length(object@C)), "and the length of grid is", as.character(length(object@grid)), "whereas the length of these lists should be equal.")
    
    errors <- c(errors, msg)
    
    msg <- NULL
  }
  
  if (length(object@B) != length(object@grid)) {
    msg <- paste("The length of B is", as.character(length(object@B)), "and the length of grid is", as.character(length(object@grid)), "whereas the length of these lists should be equal.")
    
    errors <- c(errors, msg)
    
    msg <- NULL
  }
  
  # check if the classes and values for supplied inputs are correct
  for (i in 1L:length(object@C)) {
    if (class(object@C[[i]]) != "matrix" && class(object@C[[i]]) != "numeric" && class(object@C[[i]]) != "integer" && class(object@C[[i]]) != "array") {
      msg <- paste("The class of entry", as.character(i), "of X is not matrix, numeric, or integer as is expected.")
      
      errors <- c(errors, msg)
      
      msg <- NULL
    }
  }
  
  for (i in 1L:length(object@B)) {
    if (class(object@B[[i]]) != "matrix" && class(object@B[[i]]) != "numeric" && class(object@B[[i]]) != "integer" && class(object@B[[i]]) != "list") {
      msg <- paste("The class of entry", as.character(i), "of B is not matrix, numeric, integer, or list as is expected.")
      
      errors <- c(errors, msg)
      
      msg <- NULL
    }
    
    if (class(object@grid[[i]]) == "list" && length(object@grid[[i]]) == 2) {
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) != 4) {
        msg <- paste("The length of entry", as.character(i), "of B is not four. The entry should be a basis matrix or a list of the form list(int_1,int_2,string_1,string_2) for variables taken over a two-dimensional domain where int_1,int_2>0 and string_1,string_2 = \"bspline\" or string_1,string_2 = \"fourier\".")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 4 && class(object@B[[i]][[1]]) != "integer" && class(object@B[[i]][[1]]) != "numeric") {
        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 4 && class(object@B[[i]][[1]]) == "integer" || class(object@B[[i]][[1]]) == "numeric" && object@B[[i]][[1]] <= 0) {
        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 4 && class(object@B[[i]][[2]]) != "integer" && class(object@B[[i]][[2]]) != "numeric") {
        msg <- paste("The second element in entry", as.character(i), "of B should be a positive integer.")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 4 && class(object@B[[i]][[2]]) == "integer" || class(object@B[[i]][[2]]) == "numeric" && object@B[[i]][[2]] <= 0) {
        msg <- paste("The second element in entry", as.character(i), "of B should be a positive integer.")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 4 && class(object@B[[i]][[3]]) != "character") {
        msg <- paste("The third element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 4 && class(object@B[[i]][[3]]) == "character" && object@B[[i]][[3]] != "bspline" && object@B[[i]][[3]] != "fourier") {
        msg <- paste("The third element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 4 && class(object@B[[i]][[4]]) != "character") {
        msg <- paste("The fourth element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 4 && class(object@B[[i]][[4]]) == "character" && object@B[[i]][[4]] != "bspline" && object@B[[i]][[4]] != "fourier") {
        msg <- paste("The fourth element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
    } else {
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) != 2) {
        msg <- paste("The length of entry", as.character(i), "of B is not two. The entry should be a basis matrix or a list of the form list(int,string) for variables taken over a one-dimensional domain where int>0 and string = \"bspline\" or string = \"fourier\".")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 2 && class(object@B[[i]][[1]]) != "integer" && class(object@B[[i]][[1]]) != "numeric") {
        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 2 && class(object@B[[i]][[1]]) == "integer" || class(object@B[[i]][[1]]) == "numeric" && object@B[[i]][[1]] <= 0) {
        msg <- paste("The first element in entry", as.character(i), "of B should be a positive integer.")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 2 && class(object@B[[i]][[2]]) != "character") {
        msg <- paste("The second element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      if (class(object@B[[i]]) == "list" && length(object@B[[i]]) == 2 && class(object@B[[i]][[2]]) == "character" && object@B[[i]][[2]] != "bspline" && object@B[[i]][[2]] != "fourier") {
        msg <- paste("The second element in entry", as.character(i), "of B should be a string that is either \"bspline\" or \"fourier\".")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
    }
  }
  # Check classes and supplied inputs for the grids
  for (i in 1L:length(object@grid)) {
    if (class(object@grid[[i]][1]) == "list") {
      for (k in 1:length(object@grid[[i]])) {
        if (class(object@grid[[i]][[k]]) != "matrix" && class(object@grid[[i]][[k]]) != "integer" && class(object@grid[[i]][[k]]) != "numeric") {
          msg <- paste("Element ", as.character(k), " of entry ", as.character(i), " of grid is not a matrix, numeric, or integer as is expected.", sep = "")
          errors <- c(errors, msg)
          
          msg <- NULL
        }
        
        if (class(object@grid[[i]][[k]]) == "numeric" && object@grid[[i]][[k]][1] >= object@grid[[i]][[k]][2]) {
          msg <- paste("The first element of interval ", as.character(k), " of entry ", as.character(i), " of grid should be less than the second element of that same interval.", sep = "")
          
          errors <- c(errors, msg)
          
          msg <- NULL
        }
      }
    } else {
      if (class(object@grid[[i]]) != "matrix" && class(object@grid[[i]]) != "numeric" && class(object@grid[[i]]) != "integer") {
        msg <- paste("The class of entry ", as.character(i), " of grid is not matrix, numeric, integer, or list as is expected.", sep = "")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
      
      
      if (class(object@grid[[i]]) == "numeric" && object@grid[[i]][1] >= object@grid[[i]][2]) {
        msg <- paste("The first element of entry ", as.character(i), " of grid should be less than the second element.", sep = "")
        
        errors <- c(errors, msg)
        
        msg <- NULL
      }
    }
  }
  
  if (length(errors) != 0) {
    return(errors)
  }
}




# Set the class definition
setClass("mfd", representation(C = "list", B = "list", grid = "list"), validity = check_mfd)

