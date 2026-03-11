library(testthat)

test_that("calc_size_factors computes correct size factors", {
        ## Create a mock 2x2 matrix of coverage proxies
        mock_matrix <- matrix(
                c(2, 4, 8, 16), 
                nrow = 2, 
                dimnames = list(c("tax1", "tax2"), c("sample1", "sample2"))
        )
        
        ## Calculate size factors using the median method
        sf_median <- calc_size_factors(
                x = mock_matrix, 
                method = "median", 
                verbose = FALSE
        )
        
        ## Verify the output structure and dimensions
        expect_equal(length(sf_median), 2)
        expect_true(is.numeric(sf_median))
        expect_named(sf_median, c("sample1", "sample2"))
})