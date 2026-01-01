test_that("extract_kmer_counts sums correct taxid counts", {
  kmer_string <- "123:5 456:10 123:2"
  expect_equal(basen:::extract_kmer_counts(kmer_string, 123), 7)
  expect_equal(basen:::extract_kmer_counts(kmer_string, 456), 10)
})