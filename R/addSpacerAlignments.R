
#' @inherit search_index
#' @title CHOPOFF version of crisprDesign::addSpacerAlignments
#' @param guideSet a GuideSet object from the crisprDesign package
#' @return the same guideSet object with added alignments to \code{mcols(guideSet)[["alignments"]]}
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors List
#' @export
addSpacerAlignmentsCHOPOFF <- function(guideSet, index_dir, out_file = file.path(index_dir, paste0(algorithm, "_", distance, ".csv")),
                                       algorithm = "prefixHashDB",
                                       distance = 3, validate = TRUE,
                                       chopoff_path = CHOPOFF_renviron()) {
  if (!is(guideSet, "GuideSet")) stop("Requires input of class GuideSet from the crisprDesign package",
                                      "please see the vignette for an example if you are unsure how to make one.")
  message("- Searching for alignments:")
  guide_hits <- search_index(guideSet, index_dir, out_file = out_file,
                             distance = distance, algorithm = algorithm,
                             validate = validate, chopoff_path = chopoff_path)

  message("- Summarizing alignements by distance:")
  overlap_file <- file.path(tempdir(), "overlaps.csv")
  summarize_overlaps(guide_hits, out_file = overlap_file, distance = distance,
                     validate = FALSE)
  alignments <- read.table(overlap_file, sep = ",", header = TRUE)
  colnames(alignments)[-1] <- gsub("D", "n", colnames(alignments)[-1])

  # Read in CHOPOFF hits from guide search
  message("- Inserting allignments to mcols")
  guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)

  # Convert back to crisprVerse GuideSet object
  # Make mismatch scoring object
  matching <- match(as.character(guideSet$protospacer),  alignments$guide)
  alignments_matched <- alignments[matching, ]
  mcols(guideSet) <- cbind(mcols(guideSet), alignments_matched[, -1])
  # Now make the GuideSet [["alignments]] mcol needed for scoring
  #  mcols(guideSet)[["alignments"]] # <- This is what we want to make
  matched_dups_a <- match(alignments_matched$guide, guide_hits_table$guide)
  matched_dups_g <- match(guide_hits_table$guide, alignments_matched$guide)
  ordering <- order(matched_dups_g)
  guides_hits_table_ordered <- guide_hits_table[ordering,]
  N <- table(guides_hits_table_ordered$guide)[unique(guides_hits_table_ordered$guide)]
  dups_total <- rep.int(seq(length(guideSet)), times = N)
  # Create alignment GRanges object
  alignments_grl <- GRanges(seqnames = guides_hits_table_ordered$chromosome,
                            ranges = IRanges(guides_hits_table_ordered$start, width = 1),
                            strand = guides_hits_table_ordered$strand)
  names(alignments_grl) <- names(guideSet[dups_total])
  # We should have defined spacer as: guides_hits_table_ordered$alignment_guide
  # But crisprVerse does not support this, so we use guide.

  positive_strand <- as.character(strand(alignments_grl)) == "+"
  mcols(alignments_grl) <- cbind(spacer = guides_hits_table_ordered$guide,
                                 mcols(guideSet[dups_total])[c("protospacer", "pam", "pam_site")],
                                 n_mismatches = guides_hits_table_ordered$distance,
                                 canonical = TRUE,
                                 cut_site = start(alignments_grl) + ifelse(positive_strand, -3, 3))
  alignments_grl <- split(alignments_grl, names(alignments_grl))[unique(names(alignments_grl))]
  mcols(guideSet)[["alignments"]] <- List(alignments_grl)
  return(guideSet)
}
