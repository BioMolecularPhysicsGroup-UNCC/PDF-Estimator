estimatePDF <- function(sample, pdfLength = NULL, estimationPoints = NULL,
                        lowerBound = NULL, upperBound = NULL, target = 70, 
                        lagrangeMin = 1, lagrangeMax = 200, 
                        debug = 0, outlierCutoff = 7, smooth = TRUE) {

  
  
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
  
  if (!is.null(upperBound) && !is.numeric(upperBound)) {
    stop("upperBound must be numeric")
  }
  if (!is.null(lowerBound) && !is.numeric(lowerBound)) {
    stop("lowerBound must be numeric")
  }
  if (!is.null(lowerBound) && !is.null(upperBound)) {
    if (lowerBound >= upperBound) {
      stop("lowerBound must be less than upperBound")
    }
  }

  if (!is.numeric(lagrangeMin)) {
    stop("lagrangeMin must be an integer greater than zero")
  }
  if (!is.numeric(lagrangeMax)) {
    stop("lagrangeMax must be an integer greater than zero")
  }
  if (lagrangeMin >= lagrangeMax) {
    stop("lagrangeMin must be less than or equal to lagrangeMax")
  }
 
  if (!is.numeric(target) || target < 1 || target > 100) {
    stop("target confidence percentage must be between 1 and 100")
  }
  
  if (!is.numeric(outlierCutoff) || outlierCutoff < 0) {
    stop("outlierCutoff must be a number greater than 0")
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
	
  estimationLength = vector("numeric", 1)
  isSpecifyPoints = vector("numeric", 1)
  if (is.null(estimationPoints)) {
    isSpecifyPoints[1] = FALSE
  } else {
    if(!is.vector(estimationPoints, mode = "any") || !is.numeric(estimationPoints)) {
      stop("estimationPoints must be a vector of numeric data")
    }
    isSpecifyPoints[1] = TRUE
    estimationLength[1] = length(estimationPoints)
    if (length(estimationPoints) > pdfLength) {
      pdfLength = length(estimationPoints)
    }
  }
  
  outputLength = vector("numeric", 1)
  outputLength[1] = pdfLength
  
	low = vector("numeric", 1)
	high = vector("numeric", 1)
	
	if (is.null(lowerBound)) {
		isLow = vector("numeric", 1)
		isLow[1] = FALSE
	} else {		
		low[1] = lowerBound
		isLow = vector("numeric", 1)
		isLow[1] = TRUE
	}
	if (is.null(upperBound)) {
		isHigh = vector("numeric", 1)
		isHigh[1] = FALSE
	} else {		
		high[1] = upperBound
		isHigh = vector("numeric", 1)
		isHigh[1] = TRUE
	}

	
	failedSolution = vector("numeric", 1)
	failedSolution [1] = 0
	
	threshold = vector("numeric", 1)
	threshold [1] = 0

	sqrSize = vector("numeric", 1)
	sqrSize [1] = 0
	
	distribution=.C(getNativeSymbolInfo("estimatePDF", "PDFEstimator"), 
				sample = as.double(sample), 
				inputLength = as.integer(inputLength),
				estimationPoints = as.double(estimationPoints),
				pointsPDF = as.double(vector("numeric", estimationLength)),
				estimationLength = as.integer(estimationLength),
				isSpecifyPoints = as.integer(isSpecifyPoints),
				low = as.double(low), 	
				high = as.double(high), 
				isLow = as.integer(isLow), 
				isHigh = as.integer(isHigh),
				target = as.double(target),
				lagrangeMin = as.integer(lagrangeMin),
				lagrangeMax = as.integer(lagrangeMax),
				debug = as.integer(debug),
				outlierCutoff = as.integer(outlierCutoff),
				smooth = as.integer(smooth),
				outputLength = as.integer(outputLength),
				failedSolution = as.integer(failedSolution),
				threshold = as.double(threshold),
				x = as.double(vector("numeric", outputLength)), 
				pdf = as.double(vector("numeric", outputLength)), 
				cdf = as.double(vector("numeric", outputLength)),
				sqr = as.double(vector("numeric", inputLength)),
				sqrSize = as.double(sqrSize),
				lagrange = as.double(vector("numeric", lagrangeMax)),
				r = as.double(vector("numeric", inputLength)))
	
	      class(distribution) <- "PDFe"
	      return(distribution)
}
