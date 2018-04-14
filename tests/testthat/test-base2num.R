library(Hapi)

context("base2num")

test_that("Convert A/T/C/G to 0/1", {
    ref <- rep('A',10)
    alt <- rep('G',10)
    gmtDa <- as.matrix(sample(c('A','G'), replace = TRUE, size = 50),10,5)
    gmtDa <- base2num(gmt=gmtDa, ref=ref, alt=alt)
    
    expect_equal(sum(sapply(as.numeric(unlist(gmtDa)), is.numeric)),50)
})
