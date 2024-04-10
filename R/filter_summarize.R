#' filter overlaps

#' @param distance  What is the distance for overlap filtering of off-targets. (default: 3)
#' @param detail_file Path to the file where database is stored (output of search_index function, .csv file)
#' @param out_file File path to the file where summarized output should be generated.
#' @param validate TRUE, if false, do not check that CHOPOFF path is valid
#' @param chopoff_path PATH to CHOPOFF, default install_chopoff()
#' @return path to csv output file of filtered guides
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
#'
#'filter_overlaps(distance = 3, guide_hits)
filter_overlaps <- function(distance = 3,detail_file,out_file = file.path(dirname(detail_file), paste0(sub(".csv","",basename(detail_file)), "_filter_", distance, ".csv")),
                            validate = TRUE,
                         chopoff_path = install_CHOPOFF()) {
  stopifnot(file.exists(detail_file))

  if (validate) check_exist_and_get_version(chopoff_path)

  args <- c("--distance" = distance, "--detail_file" = detail_file, "--output" = out_file)
  args <- c("filter", paste(names(args), shQuote(args)))
  system(paste(normalizePath(chopoff_path), paste(args, collapse = " ")), wait = TRUE)
  return(out_file)
}

#' summarize overlaps


#' @param detail_file Path to the file where database is stored (output of search_index function, .csv file)
#' @param out_file File path to the file where summarized output should be generated.
#' @param validate TRUE, if false, do not check that CHOPOFF path is valid
#' @param chopoff_path PATH to CHOPOFF, default install_chopoff()
#' @return path to csv output file of filtered guides
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
#'
#'summarize_overlaps(guide_hits)
summarize_overlaps <- function(distance = 3,detail_file,out_file = file.path(dirname(detail_file), paste0(sub(".csv","",basename(detail_file)), "_summarized", ".csv")),
                            validate = TRUE,
                            chopoff_path = install_CHOPOFF()) {
  stopifnot(file.exists(detail_file))

  if (validate) check_exist_and_get_version(chopoff_path)

  args <- c("--detail_file" = detail_file, "--output" = out_file)
  args <- c("summarize", paste(names(args), shQuote(args)))
  system(paste(normalizePath(chopoff_path), paste(args, collapse = " ")), wait = TRUE)
  return(out_file)
}
