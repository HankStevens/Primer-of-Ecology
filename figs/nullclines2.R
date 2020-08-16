nullclines2 <- function (deriv, xlim, ylim, parameters = NULL, system = "two.dim", 
          points = 101, col = c("blue", "cyan"), lty=2:1, add = TRUE, add.legend = TRUE, 
          state.names = c("x", "y"), ...) 
{
  if ((!is.vector(xlim)) | (length(xlim) != 2)) {
    stop("xlim is not a vector of length 2 as required")
  }
  if (xlim[2] <= xlim[1]) {
    stop("xlim[2] is less than or equal to xlim[1]")
  }
  if ((!is.vector(ylim)) | (length(ylim) != 2)) {
    stop("ylim is not a vector of length 2 as required")
  }
  if (ylim[2] <= ylim[1]) {
    stop("ylim[2] is less than or equal to ylim[1]")
  }
  if (points <= 0) {
    stop("points is less than or equal to zero")
  }
  if (!(system %in% c("one.dim", "two.dim"))) {
    stop("system must either be set to \"one.dim\" or \"two.dim\"")
  }
  if (is.vector(col) == FALSE) {
    stop("col is not a vector as required")
  }
  if (length(col) != 2) {
    if (length(col) == 1) {
      col <- rep(col, 2)
    }
    if (length(col) > 2) {
      col <- col[1:2]
    }
    print("Note: col has been reset as required")
  }
  if (!is.logical(add)) {
    stop("add must be logical")
  }
  if (!is.logical(add.legend)) {
    stop("add.legend must be logical")
  }
  x <- seq(from = xlim[1], to = xlim[2], length = points)
  y <- seq(from = ylim[1], to = ylim[2], length = points)
  dx <- matrix(0, ncol = points, nrow = points)
  dy <- matrix(0, ncol = points, nrow = points)
  if (system == "one.dim") {
    for (i in 1:points) {
      dy[1, i] <- deriv(0, setNames(c(y[i]), state.names[1]), 
                        parameters)[[1]]
    }
    for (i in 2:points) {
      dy[i, ] <- dy[1, ]
    }
    contour(x, y, dy, levels = 0, add = add, col = col[1], 
            drawlabels = FALSE, ...)
    if (add.legend == TRUE) {
      legend("bottomright", "dy/dt = 0 for all t", lty = 1, 
             lwd = 1, col = col[1])
    }
    return(list(add = add, add.legend = add.legend, col = col, 
                deriv = deriv, dy = dy, parameters = parameters, 
                points = points, system = system, x = x, xlim = xlim, 
                y = y, ylim = ylim))
  }
  else {
    for (i in 1:points) {
      for (j in 1:points) {
        df <- deriv(0, setNames(c(x[i], y[j]), state.names), 
                    parameters)
        dx[i, j] <- df[[1]][1]
        dy[i, j] <- df[[1]][2]
      }
    }
    contour(x, y, dx, levels = 0, add = add, col = col[1], lty=lty[1],
            drawlabels = FALSE, ...)
    contour(x, y, dy, levels = 0, add = TRUE, col = col[2], lty=lty[2],
            drawlabels = FALSE, ...)
    if (add.legend == TRUE) {
      legend("bottomright", c("x nullclines", "y nullclines"), 
             lty = 1, lwd = 1, col = col)
    }
    return(list(add = add, add.legend = add.legend, col = col, 
                deriv = deriv, dx = dx, dy = dy, parameters = parameters, 
                points = points, system = system, x = x, xlim = xlim, 
                y = y, ylim = ylim))
  }
}
