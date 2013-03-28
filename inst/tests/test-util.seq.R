context("util.seq")

test_that("count.aa() warns about unrecognized amino acids and performs substring operations", {
  expect_message(count.aa("ABCDEFGHIJ"), "count.aa: unrecognized amino acid code\\(s\\): B J")
  myseq <- "AAAAAGGGGG"
  expect_equal(count.aa(myseq, stop=5)[, "G"], 0)
  expect_equal(count.aa(myseq, start=6)[, "A"], 0)
  expect_equal(count.aa(myseq, start=5, stop=6)[, c("A", "G")], c(1, 1), check.attributes=FALSE)
})
