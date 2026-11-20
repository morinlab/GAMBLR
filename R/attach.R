GAMBLR.core <- c("GAMBLR.data", "GAMBLR.helpers", "GAMBLR.utils", "GAMBLR.viz", "GAMBLR.results")

GAMBLR.core_unloaded <- function() {
    search <- paste0("package:", GAMBLR.core)
    GAMBLR.core[!search %in% search()]
}

# Attach the package from the same package library it was loaded from before
GAMBLR.same_library <- function(pkg) {
    loc <- if (pkg %in% loadedNamespaces()) dirname(getNamespaceInfo(pkg, "path"))
    library(pkg, lib.loc = loc, character.only = TRUE, warn.conflicts = FALSE)
}

GAMBLR_attach <- function() {
    to_load <- GAMBLR.core_unloaded()

    suppressPackageStartupMessages(
        lapply(to_load, GAMBLR.same_library)
    )

    invisible(to_load)
}

# Path to the user's config, honouring R_CONFIG_FILE like config::get() does
GAMBLR_config_file <- function() {
    Sys.getenv("R_CONFIG_FILE", "config.yml")
}

# The template config.yml shipped with GAMBLR.results
GAMBLR_config_template <- function() {
    system.file("extdata", "config.yml", package = "GAMBLR.results")
}

# Copy the template config.yml shipped with GAMBLR.results into the project
# root. Not exported; call as GAMBLR:::use_gamblr_config(). Users should edit
# project_base and repo_base afterwards.
use_gamblr_config <- function(path = ".", overwrite = FALSE) {
    template <- GAMBLR_config_template()
    if (!nzchar(template)) {
        cli::cli_abort("The installed GAMBLR.results doesn't include a template {.file config.yml}. Update GAMBLR.results and try again.")
    }
    dest <- file.path(path, "config.yml")
    if (file.exists(dest) && !overwrite) {
        cli::cli_abort(c(
            "{.file {dest}} already exists.",
            "i" = "Use {.code overwrite = TRUE} to replace it, or compare it against {.file {template}}."
        ))
    }
    if (!file.copy(template, dest, overwrite = overwrite)) {
        cli::cli_abort("Couldn't copy {.file {template}} to {.file {dest}}.")
    }
    cli::cli_inform(c(
        "v" = "Copied the GAMBLR config template to {.file {dest}}.",
        "!" = "Edit {.field project_base} and {.field repo_base} to match your setup."
    ))
    invisible(dest)
}

# R CMD INSTALL and R CMD check load the package to test it, from a directory
# with no config.yml, so skip the config checks there
GAMBLR_installing_or_checking <- function() {
    nzchar(Sys.getenv("R_INSTALL_PKG")) || nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))
}

GAMBLR_check_config_exists <- function() {
    config_file <- GAMBLR_config_file()
    if (file.exists(config_file)) {
        return(invisible(config_file))
    }
    if (interactive() && isTRUE(utils::askYesNo(
        paste0("No ", config_file, " found in ", getwd(), ". Copy the GAMBLR template here?")
    ))) {
        return(use_gamblr_config(dirname(config_file)))
    }
    cli::cli_abort(c(
        "Can't find {.file {config_file}} in {.path {getwd()}}.",
        "i" = "GAMBLR expects {.file config.yml} in your project root. Run {.run GAMBLR:::use_gamblr_config()} to start from the template, set your working directory to the project root, or point {.envvar R_CONFIG_FILE} at your config."
    ))
}

# Compare the active config against the keys the sub-packages read, as listed
# in GAMBLR.results/inst/extdata/required_config_keys.txt
GAMBLR_missing_config_keys <- function(config_file = GAMBLR_config_file()) {
    manifest <- system.file("extdata", "required_config_keys.txt", package = "GAMBLR.results")
    if (!nzchar(manifest)) {
        return(character(0))
    }
    required <- readLines(manifest)
    config_list <- config::get(file = config_file)

    has_key <- function(key) {
        x <- config_list
        for (part in strsplit(key, "$", fixed = TRUE)[[1]]) {
            if (!is.list(x) || !part %in% names(x)) {
                return(FALSE)
            }
            x <- x[[part]]
        }
        TRUE
    }
    required[!vapply(required, has_key, logical(1))]
}

GAMBLR_check_config_keys <- function() {
    missing_keys <- GAMBLR_missing_config_keys()
    if (length(missing_keys) > 0) {
        names(missing_keys) <- rep("x", length(missing_keys))
        cli::cli_warn(c(
            "Your {.file {GAMBLR_config_file()}} (config: {.val {Sys.getenv('R_CONFIG_ACTIVE', 'default')}}) is missing {length(missing_keys)} key{?s} used by GAMBLR functions:",
            missing_keys,
            "i" = "Functions that rely on these keys will fail. Copy the missing keys from the template at {.file {GAMBLR_config_template()}}."
        ))
    }
    invisible(missing_keys)
}

GAMBLR_attach_message <- function(to_load) {
    if (length(to_load) == 0) {
        return(NULL)
    }
    
    header <- cli::rule(
        left = cli::style_bold("Welcome fellow GAMBLer! Attaching core GAMBLR packages"),
        right = paste0("GAMBLR ", package_version_h("GAMBLR"))
    )

    to_load <- sort(to_load)
    versions <- vapply(to_load, package_version_h, character(1))

    packages <- paste0(
        cli::col_green(cli::symbol$tick), " ", cli::col_yellow(format(to_load)), " ",
        cli::ansi_align(versions, max(cli::ansi_nchar(versions)))
    )

    if (length(packages) %% 2 == 1) {
        packages <- append(packages, "")
    }
    col1 <- seq_len(length(packages) / 2)
    info <- paste0(packages[col1], "     ", packages[-col1])

    paste0(header, "\n", paste(info, collapse = "\n"))
}

package_version_h <- function(pkg) {
    highlight_version(utils::packageVersion(pkg))
}

highlight_version <- function(x) {
    x <- as.character(x)

    is_dev <- function(x) {
        x <- suppressWarnings(as.numeric(x))
        !is.na(x) & x >= 9000
    }

    pieces <- strsplit(x, ".", fixed = TRUE)
    pieces <- lapply(pieces, function(x) ifelse(is_dev(x), cli::col_red(x), cli::col_blue(x)))
    vapply(pieces, paste, collapse = ".", FUN.VALUE = character(1))
}

inform_startup <- function(msg, ...) {
    if (is.null(msg)) {
        return()
    }

    rlang::inform(msg, ..., class = "packageStartupMessage")
}

.onAttach <- function(...) {
    check_config <- !GAMBLR_installing_or_checking()
    if (check_config) {
        GAMBLR_check_config_exists()
    }
    attached <- GAMBLR_attach()

    inform_startup(GAMBLR_attach_message(attached))
    logo_lines = paste(c(
      "",
      r"{  /$$$$$$     /$$$$$$    /$$      /$$   /$$$$$$$    /$$        .:::::::    }" ,
      r"{ /$$__  $$   /$$__  $$  | $$$    /$$$  | $$__  $$  | $$        .::    .::  }" ,
      r"{| $$  \__/  | $$  \ $$  | $$$$  /$$$$  | $$  \ $$  | $$        .::    .::  }" ,
      r"{| $$ /$$$$  | $$$$$$$$  | $$ $$/$$ $$  | $$$$$$$   | $$        .: .::      }" ,
      r"{| $$|_  $$  | $$__  $$  | $$  $$$| $$  | $$__  $$  | $$        .::  .::    }" ,
      r"{| $$  \ $$  | $$  | $$  | $$\  $ | $$  | $$  \ $$  | $$        .::    .::  }" ,
      r"{|  $$$$$$/  | $$  | $$  | $$ \/  | $$  | $$$$$$$/  | $$$$$$$$  .::      .::}" ,
      r"{ \______/   |__/  |__/  |__/     |__/  |_______/   |________/}",
      r"{ ~GENOMIC~~~~~~~~~~~~~OF~~~~~~~~~~~~~~~~~B-CELL~~~~~~~~~~~~~~~~~~IN~~~~~~}",
      r"{ ~~~~~~~~~~~~ANALYSIS~~~~~~MATURE~~~~~~~~~~~~~~~~~~~LYMPHOMAS~~~~~~~~~~R~}",
      ""),collpse="\n")
    packageStartupMessage(logo_lines)

    if (check_config) {
        GAMBLR_check_config_keys()
    }
}
