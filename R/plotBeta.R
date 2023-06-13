plotBeta <- function (samples, resolution = 100, 
                      xPlotRange = c(min(samples), max(samples)), 
                      sqrPlotThreshold = 2) {
 
  if (!is.numeric(samples)) {
    stop("a numeric vector of sample data is required")
  } 
  if (!is.numeric(xPlotRange)) {
    stop("xPlotRange must be a numeric range for the x-axis")
  }  
  if (!is.numeric(sqrPlotThreshold) || sqrPlotThreshold < 1) {
    stop("sqrPlotThreshold must be a number greater than or equal to 1")
  }  
  if (!is.numeric(resolution) || resolution < 50 || resolution > 500) {
    stop("resolution must be a number between 50 and 500")
  }  
  
  
  nMax = resolution
  samples = sort(unique(samples))

  n = length(samples)
  nOriginal = n
  skip = 1
  if (n > nMax) {
    skip = floor(n / nMax)
    sampleIdx = seq(1, nOriginal, by = skip)
    if (nOriginal == sampleIdx[length(sampleIdx)]) {
      samples = sort(c(samples[sampleIdx]))
    } else {
      samples = sort(c(samples[sampleIdx], samples[nOriginal]))
    }
    n = length(samples)
  }

  inc = 1 / (n - 1)
  range = seq(0, 1, inc)

  xi = length(range) - 1
  x2 = vector()
  for (i in 1:xi) {
    x2[i] = (samples[i] + samples[i + 1]) / 2
  }
  x = sort(c(x2, samples))

  zBeta = dbeta(range, 1, n)
  z = matrix(zBeta / max(zBeta))

  meanBeta = 1:n / (n + 1)
  scale = sqrt(n + 2)

  sqrRange = matrix((range - meanBeta[1]) * scale)

  for (i in 2:n) {
    zBeta = dbeta(range, i, n - i + 1)
    z = cbind(z, zBeta / max(zBeta))
    sqrRange = cbind(sqrRange, (range - meanBeta[i]) * scale)
  }

  nrows = n + n + 1
  ncols = n + n - 1
  zRotate = matrix(nrow = nrows, ncol = ncols)

  transformY = seq(min(sqrRange), max(sqrRange), length.out = ncol(zRotate) - 2)

  for (i in 1:n) {
    for (j in 1:n) {
      zRotate[i+j+1, -i+j+n] = z[i, j]
    }
  }
  for (i in 2:(nrows-1)) {
    for (j in 2:(ncols-1)) {
      if (is.na(zRotate[i, j])) {
          zRotate[i, j] = (zRotate[i-1, j] + zRotate[i+1, j]) / 2
      }
    }
  }

  zRotate = zRotate[3:nrows, 2:(ncols-1)]
  zRotate[which(is.na(zRotate))] = 0

  zlim = range(z, finite = TRUE)
  nlevels = 50
  customGrays = gray.colors(100, end = 1, rev = TRUE)
  plot(NA, xlim = xPlotRange, axes = FALSE, 
       ylim = c(-sqrPlotThreshold, sqrPlotThreshold), 
       xlab="", ylab="")
  .filled.contour(x=x, y=transformY, z=zRotate, 
                  levels=pretty(zlim, nlevels), col = customGrays)
}


