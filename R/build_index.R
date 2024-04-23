
#' Build CHOPOFF index
#'
#' Build the genome index for CHOPOFF, using either a preset (Cas9 or Cas12),
#' or user defined parameters.
#'
#' If Preset is set, some other arguments are ignored.
#' @param name How will you shortly name this database?
#' @param genome Path to the genome, either .fa or .2bit. Must have a fasta
#' index in directory, named \code{paste0(genome, ".fai")}
#' @param out_dir output directory where the database will be stored
#' @param algorithm default: prefixHashDB, alternatives: TODO: list all here!
#' @param distance Within what distance can we search for off-targets? (default: 3)
#' @param preset The preset for parameters, default "Cas9". Alternatives:
#'  Cas12 and NULL (user custom settings)
#' @param hash_length numeric, default 16.
#' @param ambig_max  How many ambiguous bases are allowed inside the guide? (default: 0)
#' @param strands c("+", "-"), search both 5' (+) and 3' (-) strands
#' @param fwd_motif  Motif of pseudo-spacer in 5'-3' that will be matched
#' on the reference (without the PAM: X).
#' For example for Cas9 it is 20*N + XXX:
#'  paste(c(rep("N", 20), rep("X", 3)), collapse = "")
#' @param fwd_pam PAM in 5'-3' that will be matched on the reference (without the X).
#' For example for Cas9 it is 20*X + NGG: paste(c(rep("X", 20), "NGG"), collapse = "")
#' @param extend3prime  Defines how off-targets will be aligned to the
#' guides and where extra nucleotides will be
#' added for alignment within distance. Whether
#' to extend in the 5' and 3' direction. Default
#' is Cas9 with extend3 = false.
#' @param validate TRUE, if false, do not check that CHOPOFF path is valid
#' @param chopoff_path PATH to CHOPOFF, default install_CHOPOFF()
#' @return invisible(NULL)
#' @export
#' @examples
#' name <- "CAS9"
#' genome <- system.file("extdata/sample_genome", "semirandom.fa", package = "crisprCHOPOFF")
#' ## Note: a fasta index ".fai" file must exist in directory of genome.
#' # You can make it with:
#' #if (!file.exists(paste0(genome, ".fai"))) {
#' #Rsamtools::indexFa(genome)
#' #}
#' out_dir_index <- file.path(tempdir(), "CHOPOFF_sample_genome")
#' build_index(name, genome, out_dir_index, distance = 2)
#'
build_index <- function(name, genome, out_dir, algorithm = "prefixHashDB",
                        distance = 3, preset = "Cas9", hash_length = 16,
                        ambig_max = 0, strands = c("+", "-"),
                        fwd_motif = "NNNNNNNNNNNNNNNNNNNNXXX",
                        fwd_pam = "XXXXXXXXXXXXXXXXXXXXNGG",
                        extend3prime = FALSE,
                        validate = TRUE,
                        chopoff_path = CHOPOFF_renviron()) {
  stopifnot(is(genome, "character") && length(genome) == 1 && file.exists(genome))
  if (!file.exists(paste0(genome, ".fai"))) {
    stop("No fasta index file found, i.e. a .fai file suffix to genome path. Please run:
         Rsamtools::indexFa(genome)")
  }

  motif <- preset
  stopifnot(length(motif) %in% c(0, 1) && (is.character(motif) || is.null(motif)))
  stopifnot(is.character(strands) && all(strands %in% c("+", "-")))
  stopifnot(is.logical(extend3prime))
  if (validate) check_exist_and_get_version(chopoff_path)
  preset_check(motif, hash_length, ambig_max, strands, extend3prime)

  args <- c("--name" = name, "--genome" = genome, "--output" = out_dir,
            "--distance" =  distance, "--motif" = motif,
            "--ambig_max" = ambig_max,
            "--fwd_motif" = fwd_motif, "--fwd_pam" = fwd_pam)

  logicals <- (c("+", "-", "FALSE") %in% c(strands, extend3prime))
  logicals <- c("--not_forward", "--not_reverse", "--extend3")[!logicals]
  algorithm <- c(algorithm, "--hash_length" = hash_length)
  algorithm <- paste(names(algorithm), algorithm)
  args <- c("build", paste(names(args), shQuote(args)), logicals, algorithm)
  system(paste(normalizePath(chopoff_path), paste(args, collapse = " ")), wait = TRUE)
  return(invisible(NULL))
}

#' Search guides in CHOPOFF index
#' @param guides character vector of guides, crisprDesign::GuideSet, or
#'  path to txt file of guides.
#'  If file, the file must be a single column, 1 guide per line of correct length,
#'  without a header or row names.
#' @param index_dir the directory of the index created during 'build_index'.
#'  Will automatically find the file, as only 1 db can exist in a folder
#'  at the same time per algorithm.
#' @param out_file Path to the file where detailed results should
#' be written.
#' @param distance  Maximum edit distance to analyze. Must be less
#' or equal to the distance that was used when
#' building db. (default: 3)
#' @param validate TRUE, if false, do not check that CHOPOFF path is valid
#' @param algorithm default: prefixHashDB, alternatives: TODO: list all here!
#' @param chopoff_path PATH to CHOPOFF, default install_chopoff()
#' @return path to csv output file of guides
#' @export
#' @examples
#' name <- "CAS9"
#' genome <- system.file("extdata/sample_genome", "semirandom.fa", package = "crisprCHOPOFF")
#' ## Note: a fasta index ".fai" file must exist in directory of genome.
#' # You can make it with:
#' #if (!file.exists(paste0(genome, ".fai"))) {
#' #Rsamtools::indexFa(genome)
#' #}
#' out_dir_index <- file.path(tempdir(), "CHOPOFF_sample_genome")
#' build_index(name, genome, out_dir_index, distance = 2, validate = FALSE)
#'
#' # Now search some guides
#' guides <- system.file("extdata/sample_genome", "guides.txt", package = "crisprCHOPOFF")
#' # Quick preview in guides:
#' guide_candidates <- read.table(guides, col.names = "guides")
#' unique(nchar(unlist(guide_candidates))) # Unique lengths of guides
#' guide_hits <- search_index(guides, out_dir_index, distance = 2, validate = FALSE)
#' guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
#' # use data.table::fread for reading in large list
#' # Subset to 0 distance hits
#' dist0 <- guide_hits_table[guide_hits_table$distance == 0,]
#' head(dist0)
#' # Which chromosomes is a specific guide found on with 0 distance hits?
#' unique(dist0[dist0$guide == "TCCGGCCTGGTTATCGAAGG",]$chromosome) # 2 chromosomes
search_index <- function(guides, index_dir, out_file = file.path(index_dir, paste0(algorithm, "_", distance, ".csv")),
                         algorithm = "prefixHashDB",
                         distance = 3, validate = TRUE,
                         chopoff_path = CHOPOFF_renviron()) {
  stopifnot(dir.exists(index_dir))
  guides <- guide_input_check(guides)
  if (validate) check_exist_and_get_version(chopoff_path)

  args <- c("--guides" = guides, "--database" = index_dir, "--output" = out_file,
            "--distance" =  distance)
  args <- c("search", paste(names(args), shQuote(args)), algorithm)
  system(paste(normalizePath(chopoff_path), paste(args, collapse = " ")), wait = TRUE)
  return(out_file)
}
