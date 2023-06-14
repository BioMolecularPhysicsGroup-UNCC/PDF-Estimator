stitchPDF <- function(sample, pdfLength = NULL, debug = 0) {
  
  sampleDim = dim(sample)
  if (!is.null(sampleDim)) {
    if (sampleDim[2] > 1) {
      print("Calling estimatePDFmv for multivariate sample", quote = FALSE)
      return(estimatePDFmv(sample, debug))
    }
  }
  
  if (is.null(sample) || !is.numeric(sample)) {
    stop("a numeric vector of sample data is required")
  }
  
  inputLength = vector("numeric", 1)
  inputLength [1] = length(sample)
  
  if (is.null(pdfLength)) {
    pdfLength = floor(200 + inputLength/200.0)
    if (pdfLength > 1500) {
      pdfLength = 1500;	
    }
  } else {
    if(!is.numeric(pdfLength)) {
      stop("pdfLength must be numeric")
    }
  }
  
  outputLength = vector("numeric", 1)
  outputLength[1] = pdfLength
	
	distribution=.C(getNativeSymbolInfo("estimatePDFstitch", "PDFEstimator"), 
				sample = as.double(sample), 
				inputLength = as.integer(inputLength),
				debug = as.integer(debug),
				outputLength = as.integer(outputLength),
				x = as.double(vector("numeric", outputLength)), 
				pdf = as.double(vector("numeric", outputLength)),
				cdf = as.double(vector("numeric", outputLength)))
	      distribution$sqrSize = inputLength
	      class(distribution) <- "PDFe"
	      return(distribution)
}
