library(crisprCHOPOFF)
library(testthat)

genome <- system.file("extdata/sample_genome", "semirandom.fa", package = "crisprCHOPOFF")
guides <- system.file("extdata/sample_genome", "guides.txt", package = "crisprCHOPOFF")
guides_cas12 <- system.file("extdata/sample_genome", "guides_cas12a.txt", package = "crisprCHOPOFF")
out_dir_index <- file.path(tempdir(), "CHOPOFF_sample_genome")
out_dir_index_cas12 <- paste0(out_dir_index, "cas12")

test_that("Build CAS9 index distances (1:2)", {
  name <- "CAS9"
  build_index(name, genome, out_dir_index, validate = FALSE, distance = 1)
  build_index(name, genome, out_dir_index, validate = FALSE, distance = 2)
  expect_in("prefixHashDB.bin", list.files(out_dir_index))

})

test_that("Build CAS12 index distances (1:2)", {
  name <- "CAS12a"
  build_index(name, genome, out_dir_index_cas12, validate = FALSE, distance = 1, preset = "Cas12a")
  build_index(name, genome, out_dir_index_cas12, validate = FALSE, distance = 2, preset = "Cas12a")
  expect_in("prefixHashDB.bin", list.files(out_dir_index_cas12))
})

test_that("Search index, all input types", {
  name <- "CAS9"
  build_index(name, genome, out_dir_index, validate = FALSE, distance = 1)
  first_guide_start_pos <- 3933045

  # File path
  guide_hits <- search_index(guides, out_dir_index, validate = FALSE, distance = 1)
  guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
  expect_equal(sum(guide_hits_table$start), first_guide_start_pos)

  # Character
  guides_vector <- unlist(read.table(guides, header = FALSE))
  guide_hits <- search_index(guides_vector, out_dir_index, validate = FALSE, distance = 1)
  guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
  expect_equal(sum(guide_hits_table$start), first_guide_start_pos)

  if (requireNamespace("crisprDesign")) {
    library(crisprDesign); library(Biostrings)
    data(SpCas9, package="crisprBase")
    seqinfo <- seqinfo(Rsamtools::FaFile(genome))
    genome(seqinfo) <- "custom"
    isCircular(seqinfo) <- rep(FALSE, length(seqinfo))
    # Make seqnames
    guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
    unique_chromosome <- guide_hits_table$distance == 0 & !duplicated(guide_hits_table$guide)
    guide_hits_table <-guide_hits_table[unique_chromosome,]
    seqnames <- guide_hits_table[match(guides_vector, guide_hits_table$guide), ]$chromosome

    guideSet_sample <-
      GuideSet(ids = paste0("guide_", seq(length(guides_vector))),
             protospacers = DNAStringSet(guides_vector),
             customSequences = guides_vector,
             seqnames = seqnames,
             targetOrigin = "customSequences",
             strand = guide_hits_table$strand,
             CrisprNuclease = SpCas9,
             seqinfo = seqinfo,
             pams = rep("NGG", length(guides_vector)),
             pam_site = guide_hits_table$start)

    guide_hits <- search_index(guideSet_sample, out_dir_index, validate = FALSE, distance = 1)
    guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
    expect_equal(sum(guide_hits_table$start), first_guide_start_pos)
  }
})



test_that("test CAS9, distance 3:", {
  name <- "CAS9"
  build_index(name, genome, out_dir_index, validate = FALSE, distance = 3)

  # Quick preview in guides:
  guide_candidates <- read.table(guides, col.names = "guides")
  unique(nchar(unlist(guide_candidates))) # Unique lengths of guides
  guide_hits <- search_index(guides, out_dir_index, validate = FALSE, distance = 3)
  guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
  expect_equal(sum(guide_hits_table$start), 9206453)
})



test_that("test CAS12a, distance 2:", {
  name <- "CAS12a"
  build_index(name, genome, out_dir_index_cas12, validate = FALSE, distance = 2, preset = "Cas12a")

  # Quick preview in guides:
  guide_candidates <- read.table(guides_cas12, col.names = "guides")
  unique(nchar(unlist(guide_candidates))) # Unique lengths of guides
  guide_hits <- search_index(guides_cas12, out_dir_index_cas12, validate = FALSE, distance = 1)

  if (file.exists(guide_hits)) {
    guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
    expect_equal(sum(guide_hits_table$start), 83314)
  }
})
