#' Check if an object is of class 'basismfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'basismfd', FALSE otherwise.
#' @seealso \code{\link{is.mvbasismfd}}, \code{\link{is.mfd}}, \code{\link{is.mvmfd}}
#' @export
is.basismfd <- function(fdobj) {
  inherits(fdobj, "basismfd")
}

#' Check if an object is of class 'mvbasismfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'mvbasismfd', FALSE otherwise.
#' @seealso \code{\link{is.basismfd}}, \code{\link{is.mfd}}, \code{\link{is.mvmfd}}
#' @export
is.mvbasismfd <- function(fdobj) {
  inherits(fdobj, "mvbasismfd")
}

#' Check if an object is of class 'mfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'mfd', FALSE otherwise.
#' @seealso \code{\link{is.mvbasismfd}}, \code{\link{is.basismfd}}, \code{\link{is.mvmfd}}
#' @export
is.mfd <- function(fdobj) {
  inherits(fdobj, "mfd")
}

#' Check if an object is of class 'mvmfd'
#'
#' @param fdobj The object to check.
#' @return TRUE if the object is of class 'mvmfd', FALSE otherwise.
#' @seealso \code{\link{is.mvbasismfd}}, \code{\link{is.mfd}}, \code{\link{is.basismfd}}
#' @export
is.mvmfd <- function(fdobj) {
  inherits(fdobj, "mvmfd")
}
