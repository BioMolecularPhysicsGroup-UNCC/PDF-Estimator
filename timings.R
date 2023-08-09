library(PDFEstimator)
library(benchden)

datas <- list(list(14,  "Matterhorn"), list(11, "Normal"), list(1, "Uniform"))
sample_sizes <- list(1000, 10000, 100000, 1000000)
times <- list()

for (pair in datas) {
    for (sample_size in sample_sizes) {
        for (i in 1:3) {
            print(pair[[1]])
            x <- rberdev(sample_size, pair[[1]])
            times <- append(times, list(pair, sample_size, system.time(
                stitchPDF(x)
            )))
        }
    }
}

print(times)