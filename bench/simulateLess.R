library("PDFEstimator")
library("benchden")
library("philentropy")

source("./getEstimate.R")
source("./compareEstimates.R")


set.seed(1)
s = c(1000000)

for (distribution in 1:28) {
    for (samples in s) {
        print(c(distribution, samples))
        write(t(compareEstimates(distribution, sampleSize = samples, trials = 1, types = c("stitch"))),
            "SimulationResults.txt", append = TRUE, ncolumns = 1)
    }
}


