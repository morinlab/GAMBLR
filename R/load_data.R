# Helper functions not for export

#' Check for a version of data to load
#'
#' This function determines if a user is requesting the latest version
#' of the bundled data or wants to access one of the earlier veresions.
#'
#' @param mode Determines which data to handle. Defaults to somatic_hypermutation_locations. Will grow with more options as more data is version tracked.
#' @param this_genome_build The genome build of the data if coordinate-based. Accepts grch37 (default) or hg38.
#'
#' @return data frame
#' @import config dplyr readr GAMBLR.data
#'
#' @examples
#' determine_version()
#' determine_version(this_genome_build = "hg38")
#'
determine_version <- function(
    mode = "somatic_hypermutation_locations",
    this_genome_build = "grch37"
){
    # Determine the latest version of the data
    # Get absolute paths to bundled files in GAMBLR.data
    all_files <- system.file(
        "extdata",
        package = "GAMBLR.data"
    ) %>%
    list.files(
        recursive = TRUE,
        full.names = TRUE
    )

    # Extract version from the path
    all_files <- gsub(
        ".*extdata/",
        "",
        all_files
    )
    all_files <- all_files[grepl(mode, all_files)]

    # Determine the highst version
    versions <- gsub(".*[/]([^/]+)[/].*", "\\1", all_files)
    versions <- versions[grep('[0-9]+', versions)]
    versions <- sort(
        numeric_version(
            versions
            )
    )

    latest_version <- max(versions)

    # Which version did the user requested in config?
    requested_version <- config::get("bundled_data_versions")[[mode]]
    # Convert to numeric value if it is a string
    if(requested_version == "_latest"){
        requested_version = latest_version
    }

    # UCSC-ize the genome build format of GAMBLR
    if(this_genome_build == "grch37"){
        this_genome_build = "GRCh37"
    }else{
        this_genome_build = "GRCh38"
    }

    # Load the user-specified version of the data
    this_object <- paste0(
        mode,
        "_",
        this_genome_build,
        "_v",
        requested_version
    )

    this_data <- eval(
        parse(
            text = paste0(
                "GAMBLR.data::",
                this_object)
            )
    )

    return(this_data)

}


grch37_ashm_regions <- determine_version(
    mode = "somatic_hypermutation_locations",
    this_genome_build = "grch37"
)

hg38_ashm_regions <- determine_version(
    mode = "somatic_hypermutation_locations",
    this_genome_build = "hg38"
)

lymphoma_genes <- GAMBLR.data::get_genes(
    entities = c("BL", "MCL", "DLBCL"),
    version = config::get("bundled_data_versions")[["lymphoma_genes"]],
    gene_format = "data.frame"
)
