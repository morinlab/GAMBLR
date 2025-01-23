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

}
