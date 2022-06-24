#' @export
is.basismfd <- function(fdobj) {
  #  check whether FDOBJ is a basismfd object
  if (inherits(fdobj, "basismfd")) return(TRUE) else return(FALSE)
}


#' @export
is.mvbasismfd <- function(fdobj) {
  #  check whether FDOBJ is a mvbasismfd object
  if (inherits(fdobj, "mvbasismfd")) return(TRUE) else return(FALSE)
}


#' @export
is.mfd <- function(fdobj) {
  #  check whether FDOBJ is a mfd object
  if (inherits(fdobj, "mfd")) return(TRUE) else return(FALSE)
}


#' @export
is.mvmfd <- function(fdobj) {
  #  check whether FDOBJ is a mvmfd object
  if (inherits(fdobj, "mvmfd")) return(TRUE) else return(FALSE)
}