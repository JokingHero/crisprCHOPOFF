
#' Build CHOPOFF index
#'
#' Build the genome index for CHOPOFF, using either a preset (Cas9 or Cas12),
#' or user defined parameters.
#' @param name How will you shortly name this database?
#' @param genome Path to the genome, either .fa or .2bit. Must have a fasta
#' index in directory, named \code{paste0(genome, ".fai")}
#' @param out_dir ""
#' @param algorithm default: prefixHashDB, alternatives: TODO: list all here!
#' @param distance ""
#' @param preset The preset for parameters, default "Cas9". Alternatives:
#'  Cas12 and NULL (user custom settings)
#' @param hash_length ""
#' @param ambig_max ""
#' @param strands c("+", "-"), search both 5' (+) and 3' (-) strands
#' @param fwd_motif ""
#' @param fwd_pam ""
#' @param extend3prime ""
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
#' build_index(name, genome, out_dir_index, validate = FALSE)
#'
build_index <- function(name, genome, out_dir, algorithm = "prefixHashDB",
                        distance = 3, preset = "Cas9", hash_length = 16,
                        ambig_max = 0, strands = c("+", "-"),
                        fwd_motif = "NNNNNNNNNNNNNNNNNNNNXXX",
                        fwd_pam = "XXXXXXXXXXXXXXXXXXXXNGG",
                        extend3prime = FALSE,
                        validate = TRUE,
                        chopoff_path = install_CHOPOFF()) {
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

#' @param guides path to txt file of guides, 1 per line of correct length
#' @param out_dir ""
#' @param distance ""
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
#' build_index(name, genome, out_dir_index, validate = FALSE)
#'
#' # Now search some guides
#' guides <- system.file("extdata/sample_genome", "guides.txt", package = "crisprCHOPOFF")
#' # Quick preview in guides:
#' guide_candidates <- read.table(guides, col.names = "guides")
#' unique(nchar(unlist(guide_candidates))) # Unique lengths of guides
#' guide_hits <- search_index(guides, out_dir_index, validate = FALSE)
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
                         chopoff_path = install_CHOPOFF()) {
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

preset_check <- function(motif, hash_length, ambig_max, strands, extend3prime) {
  if (!is.null(motif)) {
    if ((hash_length != 16) | (ambig_max != 0) | !identical(strands, c("+", "-")) | (extend3prime != FALSE)) {
      warning("User defined custom parameters have no effect if preset 'motif' is set to Cas9 or Cas12,
           set motif = NULL if you want to use custom parameters!")
    }
  }
  return(invisible(NULL))
}

#' Install backend utilities for CHOPOFF
#'
#' Will install from source:\cr
#' - Julia 1.8.5\cr
#' - CHOPOFF.jl
#' @param script the bash script to run to install the backend
#' @return The path to installed CHOPOFF.jl binary, can also be retrieved with Sys.getenv("CHOPOFF").
#' @export
install_CHOPOFF <- function(script = system.file("bash_script", "install_julia_and_CHOPOFF.sh", package = "crisprCHOPOFF")) {
  path <- Sys.getenv("CHOPOFF")
  if (path == "") {
    message("- Installing Julia 1.8.5 and CHOPOFF.jl")
    message("This will only be done once, please wait 5 minutes")
    message("If any errors occur, please see github readme or alter the script")
    message("If you are running RStudio, you need to restart RStudio after install")
    system(script, wait = TRUE)
  }
  return(path)
}
