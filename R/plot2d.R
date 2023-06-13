plot2d <- function(x, xlab = "x", ylab = "y", zlab = "PDF"){
  
  estimate = x
  
  if (estimate$variableLength != 2) {
    stop("plot2d can only be used for density data with 2 variables")
  }
  
  
  persp3D(x = estimate$x[, 1], y = estimate$x[, 2], z = estimate$pdf, 
          xlab = xlab, ylab = ylab, zlab = zlab)
  
}