library(crisprCHOPOFF)
library(testthat)

genome <- system.file("extdata/sample_genome", "semirandom.fa", package = "crisprCHOPOFF")
out_dir_index <- file.path(tempdir(), "CHOPOFF_sample_genome")

test_that("test CAS9, distance 3:" {
  name <- "CAS9"
  build_index(name, genome, out_dir_index, validate = FALSE, distance = 3)

  # Now search some guides
  guides <- system.file("extdata/sample_genome", "guides.txt", package = "crisprCHOPOFF")
  # Quick preview in guides:
  guide_candidates <- read.table(guides, col.names = "guides")
  unique(nchar(unlist(guide_candidates))) # Unique lengths of guides
  guide_hits <- search_index(guides, out_dir_index, validate = FALSE, distance = 3)
  guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
  expect_equal(sum(guide_hits_table$start), 9206453)

})
#
# test_that("test CAS12a, distance 3:" {
#   name <- "CAS12a"
#   build_index(name, genome, out_dir_index, validate = FALSE, distance = 3, preset = "Cas12a")
#
#   # Now search some guides
#   guides <- system.file("extdata/sample_genome", "guides.txt", package = "crisprCHOPOFF")
#   # Quick preview in guides:
#   guide_candidates <- read.table(guides, col.names = "guides")
#   unique(nchar(unlist(guide_candidates))) # Unique lengths of guides
#   guide_hits <- search_index(guides, out_dir_index, validate = FALSE, distance = 3)
#   guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
#   expect_equal(sum(guide_hits_table$start), 9206453)
#
# })
