
#' @inherit search_index
#' @title CHOPOFF version of crisprDesign::addSpacerAlignments
#' @param crisprDesign_message logical, default TRUE. Print the message if not input is of class GuideSet.
#' @return the same GuideSet object with added alignments to \code{mcols(guideSet)[["alignments"]]}.
#'  If input was not of class GuideSet, it returns a GRangesList.
#' @importFrom GenomicRanges GRanges strand
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors List DataFrame mcols mcols<-
#' @importFrom Biostrings DNAStringSet reverseComplement subseq width
#' @export
addSpacerAlignmentsCHOPOFF <- function(guides, index_dir, out_file = file.path(index_dir, paste0(algorithm, "_", distance, ".csv")),
                                       algorithm = "prefixHashDB",
                                       distance = 3, validate = TRUE,
                                       chopoff_path = CHOPOFF_renviron(),
                                       crisprDesign_message = TRUE) {

  message("- Searching for alignments:")
  alignment_file <- search_index(guides, index_dir, out_file = out_file,
                                 distance = distance, algorithm = algorithm,
                                 validate = validate,
                                 chopoff_path = chopoff_path)

  message("- Summarizing alignements by distance:")
  overlap_file <- file.path(tempdir(), "overlaps.csv")
  summarize_overlaps(alignment_file, out_file = overlap_file,
                     distance = distance, validate = FALSE)

  message("- Inserting allignments as GRangesList to mcols")
  if (!is(guides, "GuideSet")) {
    guides <- construct_guideSet_from_table(alignment_file, crisprDesign_message)
  }
  guides <- add_spacers_as_grl(alignment_file, overlap_file, guides)


  return(guides)
}

add_spacers_as_grl <- function(guide_file, overlap_file, guideSet = NULL) {
  # Read in CHOPOFF hits from guide search
  guide_hits_table <- read.table(guide_file, sep = ",", header = TRUE)

  alignments <- read.table(overlap_file, sep = ",", header = TRUE)
  colnames(alignments)[-1] <- gsub("D", "n", colnames(alignments)[-1])

  # Convert back to crisprVerse GuideSet object
  # Make mismatch scoring object
  construct_guideSet_from_table <- is.null(guideSet)
  if (construct_guideSet_from_table) {
    guideSet <- construct_guideSet_from_table(guide_hits_table)
  }

  matching <- match(as.character(guideSet$protospacer),  alignments$guide)
  alignments_matched <- alignments[matching, ]
  mcols(guideSet) <- cbind(mcols(guideSet), alignments_matched[, -1])

  # Now make the GuideSet [["alignments]] mcol needed for scoring
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
  alignments_grl <- List(split(alignments_grl, names(alignments_grl))[unique(names(alignments_grl))])
  mcols(guideSet)[["alignments"]] <- alignments_grl

  return(guideSet)
}

construct_guideSet_from_table <- function(guide_hits_table, crisprDesign_message = TRUE) {
  if (is.character(guide_hits_table)) {
    guide_hits_table <- read.table(guide_hits_table, sep = ",", header = TRUE)
  }
  # TODO: We can easily create correct pams, by loading biostring. Decide on This.
  if (crisprDesign_message) {
    message("We reccomend to input of class GuideSet from the crisprDesign package",
            "please see the vignette for an example if you are unsure how to make one.")
  }

  unique_guides <- guide_hits_table[!duplicated(guide_hits_table$guide),]
  guideSet <- GRanges(unique_guides$chromosome, unique_guides$start,
                      unique_guides$strand)

  protospacers <- DNAStringSet(unique_guides$guide)
  pams <- subseq(protospacers, width(protospacers)-2, width(protospacers))
  mcols(guideSet) <- DataFrame(protospacer = protospacers,
                               pam = pams,
                               pam_site = unique_guides$start)
  names(guideSet) <- paste0("spacer_", seq(nrow(unique_guides)))
  return(guideSet)
}
