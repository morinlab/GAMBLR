
# Helper functions not for export

#' Check for a version of data to load
#'
#' This function determines if a user is requesting the latest version
#' of the bundled data or wants to access one of the earlier veresions.
#'
#' @param auto_connect Set to TRUE to ensure an ssh_session is created if absent
#'
#' @return data frame
#' @import config dplyr readr git2r GAMBLR.data
#'
#' @examples
determine_version <- function(
    mode = "somatic_hypermutation_locations",
    this_genome_build = "grch37"
){
    # Determine the latest version of the data
    # Specify which repo
    repo <- "https://github.com/morinlab/GAMBLR.data"

    # Get all possible tags
    all_refs <- names(
        git2r::remote_ls(repo)
    )
    tags <- all_refs[grepl("tags", all_refs)]

    # Since they are mix of character and number they are
    # not always properly ordered. Handle tag versions here
    # to get the most recent
    tags <- gsub(
        "refs/tags/|\\^.*",
        "",
        tags
    )
    tags <- sort(
        numeric_version(
            gsub(
                "v",
                "",
                tags)
            )
        )
    latest_tag <- max(tags)

    # Did the user requested latest version?
    requested_version <- config::get("bundled_data_versions")[[mode]]

    # UCSC-ize the genome build format of GAMBLR
    if(this_genome_build == "grch37"){
        this_genome_build = "GRCh37"
    }else{
        this_genome_build = "GRCh38"
    }
    # Conditionally load the latest data or user-specified version
    if(requested_version=="_latest"){
        # Read the data from github even if user
        # did not update package in a while
        path <- paste0(
            "https://raw.githubusercontent.com/morinlab/GAMBLR.data/master/inst/extdata/",
            mode,
            "_",
            this_genome_build,
            "_v",
            latest_tag,
            ".tsv"
        )
        print(path)
        this_data <- readr::read_tsv(
            path
        )
    }else{
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

    }

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
