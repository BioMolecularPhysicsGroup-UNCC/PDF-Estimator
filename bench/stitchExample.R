sampleSize = 1000
sample = rnorm(sampleSize, 0, 1)
dist = estimatePDF(sample)
plot(dist)

stitch = stitchPDF(sample, debug = 0, pdfLength = 1000)
plot(stitch)

