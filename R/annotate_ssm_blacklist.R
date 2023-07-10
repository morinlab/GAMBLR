#' @title Annotate SSM with Blacklists
#'
#' @description Annotate and auto-drop a MAF data frame with existing blacklists.
#'
#' @details Annotate and auto-drop a MAF data frame with existing blacklists to remove variants that would be dropped during the merge process.
#' This function returns a MAF format data frame with two new columns, indicating the number of occurrences of each variant in the two blacklists.
#' Note that there are a collection of parameters to this function to improve flexibility for many applications,
#' such as `return_blacklist` (returns the used blacklist to the vector given the function, or printed to the terminal if blank).
#' For returning variants that would be dropped, one can specify `invert = TRUE`, please use with caution, this is most likely the opposite of what you want from this function.
#' Lastly, the minimum count from one of the blacklists to drop a variant is specified with `drop_threshold = 4`.
#' This function also conveniently lets you know how many variants that were dropped in the annotation process.
#'
#' @param mutations_df A data frame with mutation data.
#' @param seq_type The seq_type of your mutations if you prefer to apply only the corresponding blacklist. More than one seq_type can be specified as a vector if desired. This parameter is required.
#' @param tool_name The tool or pipeline that generated the files (should be the same for all).
#' @param tool_version The version of the tool specified under `tool_name`.
#' @param annotator_name Name of annotator, default is "vcf2maf".
#' @param annotator_version Version of annotator specified under `annotator_name`.
#' @param genome_build The genome build projection for the variants you are working with (default is grch37).
#' @param project_base Optional: A full path to the directory that your blacklist_file_pattern is relative to.
#' @param blacklist_file_template Optional: A string that contains the relative path to your blacklist file from after the project_base (i.e. results) with any wildcards surrounded with curly braces.
#' @param drop_threshold The minimum count from one of the blacklists to drop a variant.
#' @param return_blacklist Boolean parameter for returning the blacklist. Default is FALSE.
#' @param use_curated_blacklist Boolean parameter for using a curated blacklist, default is FALSE.
#' @param verbose For debugging, print out a bunch of possibly useful information.
#' @param invert USE WITH CAUTION! This returns only the variants that would be dropped in the process (opposite of what you want, probably).
#'
#' @return A MAF format data frame with two new columns indicating the number of occurrences of each variant in the two blacklists.
#'
#' @import dplyr readr tidyr
#'
#' @export
#'
#' @examples
#'
#' #annotate MAF
#' deblacklisted_maf = annotate_ssm_blacklist(grande_maf,
#'                                            seq_type = "genome",
#'                                            genome_build = "hg38")
#'
annotate_ssm_blacklist = function(mutations_df,
                                  seq_type,
                                  tool_name = "slms_3",
                                  tool_version = "1.0",
                                  annotator_name = "vcf2maf",
                                  annotator_version = "1.2",
                                  genome_build = "grch37",
                                  project_base,
                                  blacklist_file_template,
                                  drop_threshold = 4,
                                  return_blacklist = FALSE,
                                  use_curated_blacklist = FALSE,
                                  verbose = FALSE,
                                  invert = FALSE){

  if(missing(seq_type)){
    message("User must specify seq_type of the mutations to select the right blacklist file. More than one seq_type can be specified as a vector if desired.")
    return()
  }

  projection = genome_build

  if(missing(blacklist_file_template)){
    blacklist_template = check_config_value(config::get("resources")$blacklist$template)
  }else{
    blacklist_template = blacklist_file_template
  }

  if(missing(project_base)){
    project_base = check_config_value(config::get("project_base"))
  }

  if(!use_curated_blacklist){
    blacklist_files = glue::glue(blacklist_template)
    blacklist_list = list()
    for(b in blacklist_files){
      full_path = paste0(project_base,b)

      #check for missingness
      if(!file.exists(full_path)){
        print(paste("missing: ", full_path))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
      }

      lifted_blacklist = suppressMessages(readr::read_tsv(full_path, col_names = c("chrpos", "blacklist_count"), col_types = "ci"))

      lifted_blacklist = lifted_blacklist %>%
        separate(chrpos, into = c("Chromosome", "Start_Position"), sep = ":") %>%
        mutate(Start_Position = as.numeric(Start_Position))

      blacklist_list[[b]] = lifted_blacklist
    }

    combined_blacklist = do.call("rbind", blacklist_list)

    if(return_blacklist){
      return(combined_blacklist)
    }

    #join using chromosome and position
    if(verbose){
      print(head(mutations_df))
      print(head(combined_blacklist))
    }

    if(str_detect(mutations_df$Chromosome, "chr")[1]){
      combined_blacklist = mutate(combined_blacklist, Chromosome = paste0("chr", Chromosome))
    }

    mutations_df = left_join(mutations_df, combined_blacklist, by = c("Chromosome", "Start_Position")) %>%
      mutate(blacklist_count = replace_na(blacklist_count, 0))
    dropped = dplyr::filter(mutations_df, blacklist_count > drop_threshold)

    if(verbose){
      if(nrow(dropped) > 0 ){
        ndrop = length(dropped$Tumor_Sample_Barcode)
        message(paste(ndrop, "variants were dropped"))
      }else{
        message("0 variants were dropped")
      }
    }

  }else{
    repo_base = check_config_value(config::get("repo_base"))
    full_path = paste0(repo_base, check_config_value(config::get("resources")$curated_blacklist))

    additional_blacklist = glue::glue(full_path) %>%
    readr::read_tsv()

    additional_blacklist = additional_blacklist %>%
      separate(chrpos, into = c("Chromosome", "Start_Position"), sep = ":") %>%
      mutate(Start_Position = as.numeric(Start_Position))

    mutations_df = left_join(mutations_df, additional_blacklist, by = c("Chromosome", "Start_Position")) %>%
      mutate(blacklist_count = tidyr::replace_na(blacklist_count, 0))

    dropped = dplyr::filter(mutations_df, blacklist_count > drop_threshold)

    if(verbose){
      if(nrow(dropped) > 0 ){
        ndrop = length(dropped$Tumor_Sample_Barcode)
        message(paste(ndrop, "variants were dropped"))
      } else {
        message("0 variants were dropped")
      }
    }
  }

  #drop anything that exceeds our threshold but keep NA
  mutations_df = dplyr::filter(mutations_df, is.na(blacklist_count) | blacklist_count < drop_threshold)
  if(invert){
    return(dropped)
  }
  return(mutations_df)
}
