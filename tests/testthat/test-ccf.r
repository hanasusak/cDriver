context("CCF")

test_that("CCF function is handling ok inputs", {
    
    expect_that(CCF("a"), throws_error("There is no mandatory VAF column!"))
    expect_that(CCF(1), throws_error("There is no mandatory VAF column!"))
    expect_that(CCF(NA), throws_error("There is no mandatory VAF column!"))
   
    df <- data.frame(VAF=c(0.2, 'a', NA))
    expect_that(CCF(df), throws_error("VAF column is not numeric!"))
    
})


test_that("CCF gives ok output", {
    
    df <- data.frame(VAF=c(0.2, 0.7, NA))
    df2 <- data.frame(VAF=c(0.2, 0.7, NA), CCF=c(0.4, 0.7, NA))
    
    expect_that(CCF(df), equals(df2))
    expect_true("CCF" %in% colnames( CCF(sample.genes.mutect)))
})


test_that("CCF without header gives ok output", {
    
    df <- data.frame(c(0.2, 0.7, NA))
    df2 <- data.frame(VAF=c(0.2, 0.7, NA), CCF=c(0.4, 0.7, NA))
    
    expect_that(CCF(df, VAF = 1), is_equivalent_to(df2))
    expect_true("CCF" %in% colnames( CCF(sample.genes.mutect)))
})