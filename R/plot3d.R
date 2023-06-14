plot3d <- function(x, xs = c(0), ys = c(0), zs = NULL, 
                   xlab = "X1", ylab = "X2", 
                   zlab = "X3"){
  
  estimate = x
  
  if (estimate$variableLength != 3) {
    stop("plot3d can only be used for density data with 3 variables")
  }
  
  
  slice3D(x = estimate$x[, 1], y = estimate$x[, 2], z = estimate$x[, 3], 
          colvar = estimate$pdf, xs = xs, ys = ys, zs = zs, 
          xlab = xlab, ylab = ylab, zlab = zlab)
  
}