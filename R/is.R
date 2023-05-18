#' Check if an object is of class 'basismfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'basismfd', FALSE otherwise.
is.basismfd <- function(fdobj) {
  inherits(fdobj, "basismfd")
}

#' Check if an object is of class 'mvbasismfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'mvbasismfd', FALSE otherwise.
is.mvbasismfd <- function(fdobj) {
  inherits(fdobj, "mvbasismfd")
}

#' Check if an object is of class 'mfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'mfd', FALSE otherwise.
is.mfd <- function(fdobj) {
  inherits(fdobj, "mfd")
}

#' Check if an object is of class 'mvmfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'mvmfd', FALSE otherwise.
is.mvmfd <- function(fdobj) {
  inherits(fdobj, "mvmfd")
}
