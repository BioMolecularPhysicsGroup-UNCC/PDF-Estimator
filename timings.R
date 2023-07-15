library(PDFEstimator)
library(benchden)

datas <- list(list(14,  "Matterhorn"), list(11, "Normal"), list(1, "Uniform"))
times <- list()

for (pair in datas) {
    for (i in 1:3) {
        print(pair[[1]])
        x <- rberdev(10000, pair[[1]])
        times <- append(times, list(pair, system.time(
            stitchPDF(x)
        )))
    }
}