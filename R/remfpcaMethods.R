#' @export
plot.remfpca <- function(remfpca_obj, expand = 0, xlab = NULL, ylab = NULL, ...) {
  if (!(inherits(remfpca_obj, "remfpca"))) stop("Argument 'remfpca_obj' is not a remfpca object.")
  percentvar <- round(100 * cumsum(remfpca_obj$sigma) / sum(remfpca_obj$sigma), 1)
  mean_mfd <- remfpca_obj$mean_mfd
  pc_mfd <- remfpca_obj$pc_mfd
  n <- pc_mfd$nobs
  p <- pc_mfd$nvar
  if (is.null(ylab)) ylab <- paste("Variable", 1:p)
  if (is.null(xlab)) xlab <- rep("time", p)
  for (ipc in 1:n) {
    op <- par(mfrow = c(p, 1), ask = TRUE)
    on.exit(par(op))
    for (j in 1:p) {
      dimSupp <- pc_mfd[ipc, j]$basis$dimSupp
      supp <- pc_mfd[ipc, j]$basis$supp
      x_grids <- seq(supp[1, 1], supp[2, 1], len = 100)
      if (expand == 0) {
        expand <- 2 * sqrt(remfpca_obj$sigma[ipc])
      }
      width <- expand * pc_mfd[ipc, j]$eval(x_grids)
      mu <- mean_mfd[1, j]$eval(x_grids)
      ylim <- c(min(mu - expand * width, mu + expand * width), max(mu - expand * width, mu + expand * width))
      plot(x_grids, mu,
        type = "l", ylim = ylim, ylab = ylab[j], xlab = xlab[j],
        main = paste("FPC", ipc, "(", percentvar[ipc], "%)"), ...
      )
      points(x_grids, mu - expand * width, pch = "-", col = 2, ...)
      points(x_grids, mu + expand * width, pch = "+", col = 3, ...)
    }
  }
}
