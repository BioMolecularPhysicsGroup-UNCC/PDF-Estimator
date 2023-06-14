estimatePDFmv <- function(sample, debug = 0, resolution = NULL) {

  sample = data.matrix(sample)
  
  if (is.null(sample) || !is.numeric(sample)) {
    stop("a numeric vector of sample data is required")
  }

  sampleDim = dim(sample)
  
  sampleLength = vector("numeric", 1)
  sampleLength [1] = sampleDim[1] 
  
  variableLength = vector("numeric", 1)
  variableLength [1] = sampleDim[2]
  
  outputLength = vector("numeric", 1)
  
  if (is.null(resolution)) {
    resolution = ceiling((sampleLength / 100) ^ (1 / (variableLength - 1)));
    if (resolution > 100) {
      resolution = 100;
    }
    if (resolution < 2) {
      resolution = 2;
    }
  }
    
  outputLength[1] = resolution
  
	distribution=.C(getNativeSymbolInfo("estimatePDFmv", "PDFEstimator"), 
				sample = as.double(sample), 
				sampleLength = as.integer(sampleLength),
				variableLength = as.integer(variableLength),
				resolution = as.integer(resolution),
				debug = as.integer(debug),
				x = as.double(vector("numeric", outputLength * variableLength)), 
				pdf = as.double(vector("numeric", outputLength ^ variableLength)))
	
	      class(distribution) <- "mvPDFe"
	      
	      distribution$x = array(distribution$x, dim = c(resolution, variableLength))
        distribution$pdf = array(distribution$pdf, dim = rep(resolution, variableLength))
	
	      return(distribution)
}
