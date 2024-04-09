
#' Build CHOPOFF index
#' @param name How will you shortly name this database?
#' @param genome Path to the genome, either .fa or .2bit. Must have a fasta
#' index in directory, named \code{paste0(genome, ".fai")}
#' @param out_dir ""
#' @param algorithm default: prefixHashDB, alternatives: TODO: list all here!
#' @param distance ""
#' @param motif ""
#' @param hash_length ""
#' @param ambig_max ""
#' @param strands c("+", "-"), search both 5' and 3' strands
#' @param fwd_motif ""
#' @param fwd_pam ""
#' @param extend3prime ""
#' @param validate TRUE, if false, do not check that CHOPOFF path is valid
#' @param chopoff_path PATH to CHOPOFF, default "CHOPOFF"
#' @examples
#' name <- "CAS9"
#' genome <- system.file("extdata/sample_genome", "semirandom.fa", package = "crisprCHOPOFF")
#' ## Note: a fasta index ".fai" file must exist in directory of genome.
#' # You can make it with:
#' #if (!file.exists(paste0(genome, ".fai"))) {
#' #Rsamtools::indexFa(genome)
#' #}
#' out_dir <- file.path(tempdir(), "CHOPOFF_sample_genome")
#' build_index(name, genome, out_dir, validate = FALSE)
#'
build_index <- function(name, genome, out_dir, algorithm = "prefixHashDB",
                        distance = 3, motif = "Cas9", hash_length = 16,
                        ambig_max = 0, strands = c("+", "-"),
                        fwd_motif = "NNNNNNNNNNNNNNNNNNNNXXX",
                        fwd_pam = "XXXXXXXXXXXXXXXXXXXXNGG",
                        extend3prime = FALSE,
                        validate = TRUE,
                        chopoff_path = "CHOPOFF") {
  stopifnot(is.character(strands))
  stopifnot(is.logical(extend3prime))
  if (validate) check_exist_and_get_version(chopoff_path)

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
}

#' Search guides in CHOPOFF index

#' @param guides path to txt file of guides, 1 per line of correct length
#' @param out_dir ""
#' @param distance ""
#' @param validate TRUE, if false, do not check that CHOPOFF path is valid
#' @param algorithm default: prefixHashDB, alternatives: TODO: list all here!
#' @param chopoff_path PATH to CHOPOFF, default "CHOPOFF"
#' @examples
#' name <- "CAS9"
#' genome <- system.file("extdata/sample_genome", "semirandom.fa", package = "crisprCHOPOFF")
#' ## Note: a fasta index ".fai" file must exist in directory of genome.
#' # You can make it with:
#' #if (!file.exists(paste0(genome, ".fai"))) {
#' #Rsamtools::indexFa(genome)
#' #}
#' out_dir_index <- file.path(tempdir(), "CHOPOFF_sample_genome")
#' build_index(name, genome, out_dir_index, validate = FALSE)
#'
#' # Now search some guides
#' guides <- system.file("extdata/sample_genome", "guides.txt", package = "crisprCHOPOFF")
#' guide_hits <- search_index(guides, out_dir_index, validate = FALSE)
#' guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
#' # use data.table::fread for reading in large list
#' # Subset to 0 distance hits
#' dist0 <- guide_hits_table[guide_hits_table$distance == 0,]
#' dist0
#' # Which chromosomes is a specific guide found on with 0 distance hits?
#' unique(dist0[dist0$guide == "TCCGGCCTGGTTATCGAAGG",]$chromosome) # 2 chromosomes
search_index <- function(guides, index_dir, out_file = file.path(index_dir, paste0(algorithm, "_", distance, ".csv")),
                         algorithm = "prefixHashDB",
                         distance = 3, validate = TRUE,
                         chopoff_path = "~/Desktop/forks/CHOPOFF.jl/build/bin/CHOPOFF") {
  stopifnot(dir.exists(index_dir))
  if (length(guides) != 1 && is.character(guides) && file.exists(guides)) {
    stop("'guides' must be character path to single existing file!")
  }
  if (validate) check_exist_and_get_version(chopoff_path)

  args <- c("--guides" = guides, "--database" = index_dir, "--output" = out_file,
            "--distance" =  distance)
  args <- c("search", paste(names(args), shQuote(args)), algorithm)
  system(paste(normalizePath(chopoff_path), paste(args, collapse = " ")), wait = TRUE)
  return(out_file)
}

check_exist_and_get_version <- function(chopoff_path) {
  stopifnot(is.character(chopoff_path) && length(chopoff_path) == 1 && chopoff_path != "")

  # Expand if path alias and check that it exists
  path <- suppressWarnings(try(system(paste("which", chopoff_path), intern = TRUE), silent = TRUE))
  if (length(path) == 0 && attr(path, "status") == 1) {
    stop("CHOPOFF is not on path, use direct link to binary like: ./CHOPOFF.jl/build/bin/CHOPOFF")
  }

  version <- try(system(paste(chopoff_path, "--version"), intern = TRUE), silent = TRUE)

  if (is(version, "try-error")) {
    stop("Could not run existing CHOPOFF binary, try running ", path, "--version in the terminal.")
  }
  message("-- Running CHOPPOFF (version:", version,")")
  return(invisible(NULL))
}
