l <- c(1, 2, 4, 8, 11, 13, 16, 17, 19, 23, 25, 26, 27, 28, 14)

library(PDFEstimator)
library(benchden)

for (dist in l) {
    x <- rberdev(100000, dist)
    s <- stitchPDF(x)
    plot(s)
    readline()
}