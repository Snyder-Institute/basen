library(testthat)

test_that("extract_kmer_counts sums correct taxid counts", {
        ## Define a mock k-mer assignment string
        kmer_string <- "123:5 456:10 123:2"
        
        ## Test specific TaxID extraction (the user-provided test)
        expect_equal(basen:::extract_kmer_counts(kmer_string, 123), 7)
        expect_equal(basen:::extract_kmer_counts(kmer_string, 456), 10)
        
        ## Test all taxa extraction (returns a named numeric array)
        all_counts <- basen:::extract_kmer_counts(kmer_string)
        expect_equal(unname(all_counts["123"]), 7)
        expect_equal(unname(all_counts["456"]), 10)
})