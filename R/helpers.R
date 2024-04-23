#' Check valid backend
#'
#' Check that Julia + CHOPPOFF.jl is installed
#' @param chopoff_path path to CHOPOFF.jl binary
#' @return invsisible(NULL)
#' @noRd
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

#' Validate input guide
#'
#' If R object, save to tmp folder as txt file, if file validate file.
#' Does not validate the guides have correct length etc.
#' @noRd
guide_input_check <- function(guides) {
  stopifnot(length(guides) > 0)
  if (is(guides, "GuideSet")) {
    stopifnot("protospacer" %in% colnames(mcols(guides)))
    guides <- as.character(guides$protospacer)
  }
  if (is.character(guides) && length(guides) >= 1) {
    is_file_input <- length(guides) == 1 && grepl("\\.txt$", guides)
    if (!is_file_input) {
      temp_file <- paste0(tempfile(), "_guides.txt")
      write.table(guides, temp_file, row.names = FALSE, col.names = FALSE,
                  quote = FALSE)
      guides <- temp_file
    }
    if (is_file_input && !file.exists(guides)) {
      stop("File input give, but guides txt file does not exist, must be .txt extension!")
    }
  } else stop("length of guides input is 0!")
  return(guides)
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

CHOPOFF_renviron <- function(stop = TRUE) {
  path <- Sys.getenv("CHOPOFF")
  if (path == "" & stop) {
    stop("Path to Julia CHOPOFF binary was not found:",
         "If you have not installed it, run: install_CHOPOFF()
         and restart R/RStudio",
         "If you have installed it already, run in user terminal:
         echo CHOPOFF=~/bin/CHOPOFF.jl/build/bin/CHOPOFF >> ~/.Renviron",
         " where path is the your install path." )
  }
  return(path)
}

#' Install backend utilities for CHOPOFF
#'
#' Will install from source:\cr
#' - Julia 1.8.5\cr
#' - CHOPOFF.jl
#' @param script the bash script to run to install the backend
#' @param path character path, default: CHOPOFF_renviron(FALSE)
#' @return The path to installed CHOPOFF.jl binary, can also be retrieved with Sys.getenv("CHOPOFF").
#' @export
install_CHOPOFF <- function(script = system.file("bash_script", "install_julia_and_CHOPOFF.sh"),
                            path = CHOPOFF_renviron(FALSE)) {
  if (path == "") {
    message("- Installing Julia 1.8.5 and CHOPOFF.jl")
    message("This will only be done once, please wait 5 minutes")
    message("If any errors occur, please see github readme or alter the script")
    message("If you are running RStudio, you need to restart RStudio after install")
    system(script, wait = TRUE)
  }
  return(path)
}
