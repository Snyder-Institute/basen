library(testthat)

test_that("filter_low_info_reads subsets correctly based on thresholds", {
        ## Create a mock k-mer profile data.table
        mock_profile <- data.table::data.table(
                read_id = c("read1", "read2", "read3"),
                unassigned_perc = c(0.2, 0.8, 0.4),
                kmers_diversity = c(1.5, 1.0, 2.5)
        )

        ## Apply filter: unassigned <= 0.5, diversity <= 2.0
        ## Expected: read1 passes. read2 fails unassigned. read3 fails diversity.
        passed_reads <- filter_low_info_reads(
                kmer_profile = mock_profile,
                max_unassigned_perc = 0.5,
                max_shannon_diversity = 2.0
        )

        ## Verify the output explicitly matches expectations
        expect_equal(length(passed_reads), 1)
        expect_equal(passed_reads, "read1")
})