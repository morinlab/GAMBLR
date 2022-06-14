#global variables
coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")
cnames = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "NAME", "SOMATIC_SCORE", "STRAND_A", "STRAND_B", "TYPE", "FILTER", "VAF_tumour", "VAF_normal", "DP_tumour", "DP_normal", "tumour_sample_id", "normal_sample_id", "pair_status")

#' Exclude samples that have been excluded from certain analyses and drop from merges
#'
#' @param tool_name The tool or pipeline that generated the files (should be the same for all).
#'
#' @return
#' @export
#'
#' @examples
#' excluded_samp = get_excluded_samples()
#'
get_excluded_samples = function(tool_name = "slms-3"){
  base = config::get("repo_base")
  excluded_df = read_tsv(paste0(base,"config/exclude.tsv"))
  excluded_samples = dplyr::filter(excluded_df, pipeline_exclude == tool_name) %>%
    pull(sample_id)

  return(excluded_samples)
}


#' Get MAF-format data frame for more than one patient using at most one augmented_maf  per patient(i.e. unique superset of variants)
#' and combine or subset from a merged MAF (wraps get_ssm_by_samples)
#' See get_ssm_by_samples for more information
#' @param these_patient_ids A vector of sample_id that you want results for. This is the only required argument.
#' @param these_samples_metadata Optional metadata table
#' @param tool_name Only supports slms-3 currently
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38)
#' @param flavour Currently this function only supports one flavour option but this feature is meant for eventual compatibility with additional variant calling parameters/versions
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs)
#' @param basic_columns Return first 43 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a list of indexes (integer) or a list of characters (matching columns in MAF).
#' @param return_cols If set to TRUE, a vector with all available column names will be returned. Default is FALSE.
#' @param subset_from_merge Instead of merging individual MAFs, the data will be subset from a pre-merged MAF of samples with the specified seq_type
#'
#' @return
#' @export
#'
#' @examples
#' patients = c("00-14595", "00-15201", "01-12047")
#' patients_ssm = get_ssm_by_patients(these_patient_ids = patients, seq_type = "genome", min_read_support = 3, basic_columns = TRUE, subset_from_merge = FALSE)
#'
get_ssm_by_patients = function(these_patient_ids,
                               these_samples_metadata,
                               tool_name = "slms-3",
                               projection = "grch37",
                               seq_type = "genome",
                               flavour = "clustered",
                               min_read_support = 3,
                               basic_columns = TRUE,
                               maf_cols = NULL,
                               return_cols = FALSE,
                               subset_from_merge = TRUE,
                               augmented = TRUE,
                               ssh_session){
  if(!subset_from_merge){
    message("WARNING: on-the-fly merges can be extremely slow and consume a lot of memory. Use at your own risk. ")
  }
  augmented = TRUE
  #always requires augmented MAFs to ensure all variants from the patient are included
  to_exclude = get_excluded_samples(tool_name)
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata() %>%
      dplyr::filter(patient_id %in% these_patient_ids) %>%
      dplyr::filter(!sample_id %in% to_exclude)
  }else{
    these_samples_metadata = these_samples_metadata %>%
      dplyr::filter(patient_id %in% these_patient_ids) %>%
      dplyr::filter(!sample_id %in% to_exclude)
  }
  # Removed because this drops all but one sample for the patient!
  #these_samples_metadata = group_by(these_samples_metadata, patient_id) %>%
  #  slice_head() %>%
  #  ungroup()

  these_sample_ids = pull(these_samples_metadata, sample_id)

  return(get_ssm_by_samples(these_sample_ids = these_sample_ids,
                            these_samples_metadata = these_samples_metadata,
                            tool_name = tool_name,
                            projection = projection,
                            seq_type = seq_type,
                            flavour = flavour,
                            min_read_support = min_read_support,
                            subset_from_merge = subset_from_merge,
                            augmented = augmented,
                            ssh_session = ssh_session,
                            basic_columns = basic_columns,
                            maf_cols = maf_cols,
                            return_cols = return_cols))
}


#' Get MAF-format data frame for more than one sample and combine together (wraps get_ssm_by_sample)
#' See get_ssm_by_sample for more information
#' @param these_sample_ids A vector of sample_id that you want results for. This is the only required argument.
#' @param these_samples_metadata Optional metadata table
#' @param tool_name Only supports slms-3 currently
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38)
#' @param flavour Currently this function only supports one flavour option but this feature is meant for eventual compatability with additional variant calling parameters/versions
#' @param these_genes A vector of genes to subset ssm to.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs)
#' @param basic_columns Return first 45 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the suer can specify what columns to be returned within the MAF. This parameter can either be a list of indexes (integer) or a list of characters.
#' @param return_cols If set to TRUE, a vector with all avaialble column names will be returned. Default is FALSE.
#' @param subset_from_merge Instead of merging individual MAFs, the data will be subset from a pre-merged MAF of samples with the specified seq_type
#' @param BETA optional argument to supply active ssh session connection for remote transfers
#'
#' @return data frame in MAF format.
#' @export
#'
#' @examples
#' merged_maf_force_unmatched = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"))
#'
get_ssm_by_samples = function(these_sample_ids,
                              these_samples_metadata,
                              tool_name = "slms-3",
                              projection = "grch37",
                              seq_type = "genome",
                              flavour = "clustered",
                              these_genes,
                              min_read_support = 3,
                              basic_columns = TRUE,
                              maf_cols = NULL,
                              return_cols = FALSE,
                              subset_from_merge = TRUE,
                              augmented = TRUE,
                              ssh_session){
  if(!subset_from_merge){
    message("WARNING: on-the-fly merges can be extremely slow and consume a lot of memory. Use at your own risk. ")
  }
  to_exclude = get_excluded_samples(tool_name)

  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata() %>%
      dplyr::filter(sample_id %in% these_sample_ids) %>%
      dplyr::filter(!sample_id %in% to_exclude)
  }else{
    these_samples_metadata = these_samples_metadata %>%
      dplyr::filter(sample_id %in% these_sample_ids) %>%
      dplyr::filter(!sample_id %in% to_exclude)
  }
  #ensure we only have sample_id that are in the remaining metadata (no excluded/unavailable samples)
  these_sample_ids = these_sample_ids[which(these_sample_ids %in% these_samples_metadata$sample_id)]

  if(flavour=="legacy"){
    warning("I lied. Access to the old variant calls is not currently supported in this function")
    # TODO: implement loading of the old merged MAF under icgc_dart... vcf2maf-1.2 ..level_3 as per the other from_flatfile functions
    return()

  }else if(flavour=="clustered"){
    if(subset_from_merge && !augmented){
      maf_template = config::get("results_flatfiles")$ssm$template$merged$deblacklisted
      maf_path = glue::glue(maf_template)
      full_maf_path =  paste0(config::get("project_base"), maf_path)
      message(paste("using existing merge:", full_maf_path))
      maf_df_merge = fread_maf(full_maf_path) %>%
        dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids)
      #subset maf to only include first 43 columns (default)
      if(basic_columns){maf_df_merge = dplyr::select(maf_df_merge, c(1:45))}
      #subset maf to a specific set of columns (defined in maf_cols)
      if(!is.null(maf_cols) && !basic_columns){maf_df_merge = dplyr::select(maf_df_merge, all_of(maf_cols))}
      #print all available columns
      if(!basic_columns && return_cols){print(colnames(maf_df_merge))}
    }

    if(subset_from_merge && augmented){
      maf_template = config::get("results_flatfiles")$ssm$template$merged$augmented
      maf_path = glue::glue(maf_template)
      full_maf_path =  paste0(config::get("project_base"), maf_path)
      message(paste("using existing merge:", full_maf_path))
      maf_df_merge = fread_maf(full_maf_path) %>%
        dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
        dplyr::filter(t_alt_count >= min_read_support)
      #subset maf to only include first 43 columns (default)
      if(basic_columns){maf_df_merge = dplyr::select(maf_df_merge, c(1:45))}
      #subset maf to a specific set of columns (defined in maf_cols)
      if(!is.null(maf_cols) && !basic_columns){maf_df_merge = dplyr::select(maf_df_merge, all_of(maf_cols))}
      #print all available columns
      if(!basic_columns && return_cols){print(colnames(maf_df_merge))}
    }

    if(!subset_from_merge){
      maf_df_list = list()
      for(this_sample in these_sample_ids){
        maf_df = get_ssm_by_sample(
          this_sample_id = this_sample,
          these_samples_metadata = these_samples_metadata,
          tool_name = tool_name,
          projection = projection,
          augmented = augmented,
          flavour = flavour,
          min_read_support = min_read_support,
          basic_columns = basic_columns,
          maf_cols = maf_cols,
          return_cols = return_cols,
          verbose = FALSE,
          ssh_session = ssh_session
        )
        maf_df_list[[this_sample]]=maf_df
      }
      maf_df_merge = bind_rows(maf_df_list)
    }
  }

  if(!missing(these_genes)){
    maf_df_merge = maf_df_merge %>%
      dplyr::filter(Hugo_Symbol %in% these_genes)
  }
  return(maf_df_merge)
}



#' Get the ssms (i.e. load MAF) for a single sample. This was implemented to allow flexibility because
#' there are some samples that we may want to use a different set of variants than those in the main GAMBL merge.
#' The current use case is to allow a force_unmatched output to be used to replace the SSMs from the merge for samples
#' with known contamination in the normal. This will also be useful to apply a blacklist to individual MAFs when coupled with
#' annotate_ssm_blacklist.
#'
#' @param this_sample_id Required. The sample_id you want the data from.
#' @param these_samples_metadata Either a single row or entire metadata table containing your sample_id
#' @param tool_name The name of the variant calling pipeline (currently only slms-3 is supported)
#' @param projection The projection genome build. Supports hg38 and grch37.
#' @param these_genes A vector of genes to subset ssm to.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF
#' @param flavour Currently this function only supports one flavour option but this feature is meant for eventual compatability with additional variant calling parameters/versions
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs)
#' @param basic_columns Return first 43 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the suer can specify what columns to be returned within the MAF. This parameter can either be a list of indexes (integer) or a list of characters.
#' @param return_cols If set to TRUE, a vector with all avaialble column names will be returned. Default is FALSE.
#' @param verbose Enable for debugging/noisier output
#' @param ssh_session BETA feature! pass active ssh session object.
#' If specified, the function will assume the user is not on the network and will temporarily copy the file locally.
#'
#' @return data frame in MAF format.
#' @export
#'
#' @examples
#' ssm_sample = get_ssm_by_sample(this_sample_id = "HTMCP-01-06-00485-01A-01D", tool_name = "slims-3", projection = "grch37")
#'
get_ssm_by_sample = function(this_sample_id,
                             these_samples_metadata,
                             tool_name = "slms-3",
                             projection = "grch37",
                             these_genes,
                             augmented = TRUE,
                             flavour = "clustered",
                             min_read_support = 3,
                             basic_columns = TRUE,
                             maf_cols = NULL,
                             return_cols = FALSE,
                             verbose = FALSE,
                             verbose = FALSE,
                             ssh_session){

  #figure out which unix_group this sample belongs to
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata() %>%
      dplyr::filter(sample_id == this_sample_id)

  }else{
    these_samples_metadata = these_samples_metadata %>%
      dplyr::filter(sample_id == this_sample_id)
  }
  sample_id = this_sample_id
  tumour_sample_id = sample_id
  unix_group = pull(these_samples_metadata, unix_group)
  genome_build = pull(these_samples_metadata, genome_build)
  target_builds = projection
  seq_type = pull(these_samples_metadata, seq_type)
  pair_status = pull(these_samples_metadata, pairing_status)
  if(verbose){
    print(paste("group:", unix_group, "genome:", genome_build))
  }
  # Get unmatched normal if necessary. This is done using the unmatched normals that were added to the GAMBLR config.
  # That will need to be kept up to date if/when any new references are added.
  if(pair_status == "unmatched"){
    normal_sample_id = config::get("unmatched_normal_ids")[[unix_group]][[seq_type]][[genome_build]]

  }else{
    normal_sample_id = pull(these_samples_metadata, normal_sample_id)
  }
  base_path = ""
  if(flavour == "legacy"){
    warning("Access to the old variant calls is not currently supported in this function")
    warning("Use get_ssm_by_samples to access the legacy flavour")
    # To be fixed maybe if we decide it's needed.
    # Implementation of backwards compatability will be a lot harder because of the old naming scheme.
    return()
  }else if(flavour == "clustered"){
    vcf_base_name = "slms-3.final"
    path_template = config::get("results_flatfiles",config="default")$ssm$template$clustered$deblacklisted
    path_complete = unname(unlist(glue::glue(path_template)))
    full_maf_path = paste0(config::get("project_base",config="default"), path_complete)
    local_full_maf_path = paste0(config::get("project_base"), path_complete)
    if(augmented){
      path_template = config::get("results_flatfiles",config="default")$ssm$template$clustered$augmented
      path_complete = unname(unlist(glue::glue(path_template)))
      aug_maf_path = paste0(config::get("project_base",config="default"), path_complete)
      local_aug_maf_path = paste0(config::get("project_base"), path_complete)
    }
  }else{
    warning("Currently the only flavour available to this function is 'clustered'")
  }
  if(!missing(ssh_session)){
    #check if file exists
    status = ssh_exec_internal(ssh_session,command=paste("stat",aug_maf_path),error=F)$status
    aug_maf_path = paste0(aug_maf_path,".gz")
    local_aug_maf_path = paste0(local_aug_maf_path,".gz")
    full_maf_path = paste0(full_maf_path,".gz")
    local_full_maf_path = paste0(local_full_maf_path,".gz")
    # first check if we already have a local copy
    # Load data from local copy or get a local copy from the remote path first
    if(status==0){
      if(verbose){
        message("found:",aug_maf_path)
        message("local home:",local_aug_maf_path)
      }
      dirN = dirname(local_aug_maf_path)

      suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
      if(!file.exists(local_aug_maf_path)){

        scp_download(ssh_session,aug_maf_path,dirN)
      }
      sample_ssm = fread_maf(local_aug_maf_path) %>%
      dplyr::filter(t_alt_count >= min_read_support)
    }else{
      if(verbose){
        message("will use",full_maf_path)
        message("local home:",full_maf_path,local_full_maf_path)
      }
      dirN = dirname(local_full_maf_path)

      suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
      if(!file.exists(local_full_maf_path)){

        scp_download(ssh_session,full_maf_path,dirN)
      }
      sample_ssm = fread_maf(local_full_maf_path)
    }
  }else if(augmented && file.exists(aug_maf_path)){
    full_maf_path = aug_maf_path
    sample_ssm = fread_maf(full_maf_path)
    if(min_read_support){
      # drop poorly supported reads but only from augmented MAF
      sample_ssm = dplyr::filter(sample_ssm, t_alt_count >= min_read_support)
    }
  }else{
    if(!file.exists(full_maf_path)){
      message(paste("ERROR: file does not exist", full_maf_path))
      return()
    }
    sample_ssm = fread_maf(full_maf_path)
  }

  if(!missing(these_genes)){
    sample_ssm = sample_ssm %>%
      dplyr::filter(Hugo_Symbol %in% these_genes)
  }

  #subset maf to only include first 43 columns (default)
  if(basic_columns){
    sample_ssm = dplyr::select(sample_ssm, c(1:45))
    }

  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    sample_ssm = dplyr::select(sample_ssm, all_of(maf_cols))
    }

  #print all available columns
  if(!basic_columns && return_cols){
    print(colnames(sample_ssm))
    }

  return(sample_ssm)
}


#' Helper function to find the production merge for a pipeline and restrict to the right file based on file permissions.
#'
#' @param tool_name Lowercase name of the tool (e.g. manta, slms-3).
#' @param projection Which genome you want your results projected (lifted) to. Currently only grch37 is supported and is default.
#' @param seq_type The seq_type you want back (currently only genome is supported/available).
#'
#' @return String that holds filepath, tool name, sequence type, projection and extension.
#' @export
#'
#' @examples
#' merged = get_merged_result("manta", "grch37", "genome")
#'
get_merged_result = function(tool_name,
                             projection = "grch37",
                             seq_type = "genome"){

  if(projection != "grch37"){
  message("Currently, only grch37 is supported")
  return()
  }

  base_path = config::get("project_base")
  gambl_only = paste0(base_path, "gambl/gamblr/02-merge/", tool_name, "/", seq_type, "/")
  gambl_plus = paste0(base_path, "all_the_things/", tool_name, "/", seq_type, "--gambl,icgc_dart/")
  if(tool_name == "manta"){
    extension = ".bedpe"
  }
  gambl_only = paste0(gambl_only, "all_", tool_name, "_merged_", projection, extension)
  gambl_plus = paste0(gambl_plus, "all_", tool_name, "_merged_", projection, extension)
  permissions = file.access(gambl_plus, 4)
  if(permissions == -1 ){
    message("restricting to non-ICGC data")
    return(gambl_only)
  }else{
    return(gambl_plus)
  }
}


#' Get GAMBL metadata.
#'
#' @param seq_type_filter Filtering criteria (default: all genomes)
#' @param tissue_status_filter Filtering criteria (default: only tumour genomes, can be "mrna" or "any" for the superset of cases)
#' @param case_set optional short name for a pre-defined set of cases avoiding any
#' @param remove_benchmarking By default the FFPE benchmarking duplicate samples will be dropped
#' @param sample_flatfile Optionally provide the full path to a samples table to use instead of the default
#' @param biopsy_flatfile Optionally provide the full path to a biopsy table to use instead of the default
#' @param with_outcomes Optionally join to gambl outcome data
#' @param only_available If TRUE, will remove samples with FALSE or NA in the bam_available column (default: TRUE)
#' @param from_flatfile New default is to use the metadata in the flatfiles from your clone of the repo. Can be over-ridden to use the database
#' @param seq_type_priority For duplicate sample_id with different seq_type available, the metadata will prioritize this seq_type and drop the others
#'
#' embargoed cases (current options: 'BLGSP-study', 'FL-study', 'DLBCL-study', 'FL-DLBCL-study', 'FL-DLBCL-all', 'DLBCL-unembargoed', 'BL-DLBCL-manuscript', 'MCL','MCL-CLL')
#'
#' @return A data frame with metadata for each biopsy in GAMBL
#' @export
#' @import tidyverse DBI RMariaDB dbplyr data.table
#'
#' @examples
#' basic usage
#' my_metadata = get_gambl_metadata()
#' use pre-defined custom sample sets
#' only_blgsp_metadata = get_gambl_metadata(case_set="BLGSP-study")
#' override default filters and request metadata for samples other than tumour genomes, e.g. also get the normals
#' only_normal_metadata = get_gambl_metadata(tissue_status_filter = c('tumour','normal'))
#' non_duplicated_genome_and_capture = get_gambl_metadata(seq_type_filter=c('genome','capture'),seq_type_priority="genome")
get_gambl_metadata = function(seq_type_filter = "genome",
                              tissue_status_filter = c("tumour"),
                              case_set,
                              remove_benchmarking = TRUE,
                              with_outcomes = TRUE,
                              from_flatfile = TRUE,
                              sample_flatfile = "",
                              biopsy_flatfile = "",
                              only_available = TRUE,
                              seq_type_priority="genome"){

  outcome_table = get_gambl_outcomes(from_flatfile = from_flatfile)

  if(from_flatfile){
    base = config::get("repo_base")
    if(sample_flatfile == ""){
      sample_flatfile = paste0(base, config::get("table_flatfiles")$samples)
    }
    if(biopsy_flatfile==""){
      biopsy_flatfile = paste0(base, config::get("table_flatfiles")$biopsies)
    }
    sample_meta = suppressMessages(read_tsv(sample_flatfile, guess_max = 100000))
    biopsy_meta = suppressMessages(read_tsv(biopsy_flatfile, guess_max = 100000))

  }else{
    db = config::get("database_name")
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
    sample_meta = dplyr::tbl(con, "sample_metadata") %>%
      as.data.frame()

    biopsy_meta = dplyr::tbl(con, "biopsy_metadata") %>%
      as.data.frame()

    DBI::dbDisconnect(con)
  }

  # Conditionally remove samples without bam_available == TRUE
  if(only_available == TRUE){
    sample_meta = dplyr::filter(sample_meta, bam_available %in% c(1, "TRUE"))
  }
  sample_meta_normal_genomes =  sample_meta %>%
    dplyr::filter(seq_type %in% seq_type_filter & tissue_status == "normal") %>%
    dplyr::select(patient_id, sample_id, seq_type, genome_build) %>% as.data.frame() %>%
    #dplyr::select(patient_id, sample_id, seq_type) %>% as.data.frame() %>%
    dplyr::rename("normal_sample_id" = "sample_id")
#print(head(sample_meta_normal_genomes))
  sample_meta = sample_meta %>%
    dplyr::filter(seq_type %in% seq_type_filter & tissue_status %in% tissue_status_filter & bam_available %in% c(1,"TRUE")) %>%
    dplyr::select(-sex)


  #if we only care about genomes, we can drop/filter anything that isn't a tumour genome
  #The key for joining this table to the mutation information is to use sample_id. Think of this as equivalent to a library_id. It will differ depending on what assay was done to the sample.
  biopsy_meta = biopsy_meta %>%
    dplyr::select(-patient_id) %>%
    dplyr::select(-pathology) %>%
    dplyr::select(-time_point) %>%
    dplyr::select(-EBV_status_inf) #drop duplicated columns

  all_meta = dplyr::left_join(sample_meta, biopsy_meta, by = "biopsy_id") %>%
    as.data.frame()

  all_meta = all_meta %>%
    mutate(bcl2_ba = ifelse(bcl2_ba == "POS_BCC", "POS", bcl2_ba))

  if(!"mrna" %in% seq_type_filter & length(tissue_status_filter) == 1 & tissue_status_filter[1] == "tumour"){
    #join back the matched normal genome
    #all_meta = left_join(all_meta, sample_meta_normal_genomes, by=c("patient_id", "seq_type"))
    all_meta = left_join(all_meta, sample_meta_normal_genomes, by=c("patient_id", "seq_type","genome_build"))
    all_meta = all_meta %>%
      mutate(pairing_status = case_when(is.na(normal_sample_id)~"unmatched", TRUE~"matched"))
  }
  #all_meta[all_meta$pathology == "B-cell unclassified","pathology"] = "HGBL"  #TODO fix this in the metadata
  if(remove_benchmarking){
    all_meta = all_meta %>%
      dplyr::filter(cohort != "FFPE_Benchmarking")
  }
  if("any" %in% seq_type_filter){
    #this option may not work. To be deprecated, probably

   all_meta = all_meta %>%
    arrange(seq_type) %>%
    group_by(biopsy_id) %>%
    slice_head()
  }
  all_meta = add_icgc_metadata(all_meta) %>%
    mutate(consensus_pathology = case_when(ICGC_PATH == "FL-DLBCL" ~ "COM",
                                           ICGC_PATH == "DH-BL" ~ pathology,
                                           ICGC_PATH == "FL" | ICGC_PATH == "DLBCL" ~ ICGC_PATH,
                                           pathology == "COMFL" ~ "COM",
                                           TRUE ~ pathology))

  all_meta = unique(all_meta) #something in the ICGC code is causing this. Need to figure out what #should this be posted as an issue on Github?
  if(!missing(case_set)){
    # This functionality is meant to eventually replace the hard-coded case sets
    case_set_path = config::get("sample_sets")$default
    full_case_set_path =  paste0(config::get("repo_base"), case_set_path)
    if (file.exists(full_case_set_path)) {
      full_case_set = suppressMessages(read_tsv(full_case_set_path))
    } else {
      message(paste("Warning: case set is requested, but the case set file", full_case_set_path, "is not found."))
      message("Defaulting to pre-defined case sets")
    }

    # pre-defined case sets
    if(case_set == "MCL"){
      all_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("MCL"))

    }else if(case_set == "MCL-CLL"){
      all_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("MCL", "CLL")) %>%
        dplyr::filter(cohort != "CLL_LSARP_Trios")
    }else if(case_set == "tFL-study"){
      #update all DLBCLs in this file to indicate they're transformations
      transformed_manual = read_tsv("/projects/rmorin/projects/gambl-repos/gambl-rmorin/data/metadata/raw_metadata/gambl_tFL_manual.tsv")

      all_meta = left_join(all_meta, transformed_manual)
      fl_meta_kridel = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        group_by(patient_id) %>%
        mutate(FL = sum(pathology == "FL"), DLBCL = sum(pathology %in% c("COM", "DLBCL", "COMFL"))) %>%
        mutate(transformed = ifelse(FL > 0 & DLBCL > 0, TRUE, FALSE))  %>%
        mutate(analysis_cohort = case_when(consensus_pathology == "FL" & transformed == TRUE ~ "pre-HT",
                                           consensus_pathology == "DLBCL" & transformed == TRUE ~ "ignore",
                                           TRUE ~ "no-HT")) %>%
        dplyr::filter(cohort == "FL_Kridel") %>%
        dplyr::filter((analysis_cohort == "no-HT" & time_point == "A")|(analysis_cohort == "pre-HT")) %>%
        dplyr::select(-transformed, -FL, -DLBCL)

      dlbcl_meta_kridel = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        group_by(patient_id) %>%
        mutate(FL = sum(pathology == "FL"), DLBCL = sum(pathology %in% c("COM", "DLBCL", "COMFL"))) %>%
        mutate(transformed = ifelse(FL > 0 & DLBCL > 0, TRUE, FALSE))  %>%
        mutate(analysis_cohort = case_when(consensus_pathology == "FL" & transformed == TRUE ~ "post-HT",
                                           consensus_pathology == "DLBCL" & transformed == TRUE ~ "ignore",
                                           TRUE ~ "post-HT")) %>%
        dplyr::filter((analysis_cohort == "post-HT" & time_point == "B"))

      fl_meta_other = all_meta %>% dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        dplyr::filter(cohort != "FL_Kridel") %>%
        dplyr::filter((consensus_pathology %in% c("FL", "COM"))) %>% mutate(analysis_cohort = consensus_pathology)

      gambl_transformations = suppressMessages(read_delim("/projects/rmorin/projects/gambl-repos/gambl-rmorin/data/metadata/raw_metadata/gambl_transformation.txt", delim = " ")) %>%
        dplyr::filter(code_transf == 1) %>%
        group_by(res_id) %>%
        slice_head()

      fl_transformation_meta = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-rmorin/shared/gambl_fl_transformed.tsv"))
      transformed_cases = pull(gambl_transformations, res_id)
      fl_meta_other[which(fl_meta_other$patient_id %in% transformed_cases), "analysis_cohort"] = "pre-HT"
      fl_meta_other = mutate(fl_meta_other, analysis_cohort = ifelse(analysis_cohort == "FL", "no-HT", analysis_cohort))

      #Finally over-ride analysis cohort with the outcome of clinical review
      dlbcl_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP", "FL_Kridel", "FFPE_Benchmarking")) %>%
        dplyr::filter(consensus_pathology == "DLBCL") %>%
        mutate(analysis_cohort = "denovo-DLBCL")

      all_meta  = bind_rows(dlbcl_meta, dlbcl_meta_kridel, fl_meta_kridel, fl_meta_other) %>%
        unique()

      curated = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-rmorin/data/metadata/raw_metadata/clin_review_fl.tsv"))

      all_meta = left_join(all_meta, curated) %>%
        mutate(analysis_cohort = ifelse(is.na(clinical_review), analysis_cohort, clinical_review))

      all_meta[which(all_meta$is_tFL == 1), "analysis_cohort"] = "post-HT"
    } else if(case_set == "FL-DLBCL-study"){

      #get FL cases and DLBCL cases not in special/embargoed cohorts
      fl_meta_kridel = all_meta %>% dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        group_by(patient_id) %>%
        mutate(FL = sum(pathology == "FL"), DLBCL = sum(pathology %in% c("COM", "DLBCL", "COMFL"))) %>%
        mutate(transformed = ifelse(FL > 0 & DLBCL > 0, TRUE, FALSE)) %>%
        mutate(analysis_cohort=case_when(consensus_pathology == "FL" & transformed == TRUE ~ "tFL", consensus_pathology == "DLBCL" & transformed == TRUE ~ "ignore", TRUE ~ "FL")) %>%
        dplyr::filter(cohort == "FL_Kridel") %>%
        dplyr::filter((analysis_cohort == "FL" & time_point == "A")|(analysis_cohort == "tFL")) %>%
        dplyr::select(-transformed, -FL, -DLBCL)

      fl_meta_other = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP")) %>%
        dplyr::filter(cohort != "FL_Kridel") %>%
        dplyr::filter((consensus_pathology %in% c("FL", "COM"))) %>% mutate(analysis_cohort = consensus_pathology)

      fl_transformation_meta = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-rmorin/shared/gambl_fl_transformed.tsv"))
      transformed_cases = fl_transformation_meta %>%
        dplyr::filter(!is.na(PATHa.tr)) %>%
        pull(patient_id)

      fl_meta_other[which(fl_meta_other$patient_id %in% transformed_cases), "analysis_cohort"] = "tFL"

      dlbcl_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL", "COM")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP", "FL_Kridel", "FFPE_Benchmarking")) %>%
        dplyr::filter(consensus_pathology == "DLBCL" & COO_final == "GCB") %>%
        mutate(analysis_cohort = "DLBCL")

      all_meta = bind_rows(dlbcl_meta, fl_meta_kridel, fl_meta_other) %>%
      unique()

    }else if(case_set == "FL-study"){

      #get FL cases and DLBCL cases not in special/embargoed cohorts
      all_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios")) %>%
        group_by(patient_id) %>%
        arrange(patient_id, pathology)  %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::filter(pathology == "FL")

    }else if(case_set == "DLBCL-study"){

      #get FL cases and DLBCL cases not in special/embargoed cohorts
      all_meta = all_meta %>%
        dplyr::filter(consensus_pathology %in% c("FL", "DLBCL")) %>%
        dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios")) %>%
        group_by(patient_id) %>%
        arrange(patient_id, pathology)  %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::filter(consensus_pathology %in% c("DLBCL", "FL"))

    }else if(case_set == "DLBCL-unembargoed"){
      #get DLBCL cases not in special/embargoed cohorts
      all_meta = all_meta %>%
      dplyr::filter(consensus_pathology %in% c("DLBCL", "COM")) %>%
      dplyr::filter(!cohort %in% c("DLBCL_ctDNA", "DLBCL_BLGSP", "LLMPP_P01", "DLBCL_LSARP_Trios", "DLBCL_HTMCP"))

    }else if(case_set == "BLGSP-study"){

      #get BL cases minus duplicates (i.e. drop benchmarking cases)
      all_meta = all_meta %>%
        dplyr::filter(cohort %in% c("BL_Adult", "BL_cell_lines", "BL_ICGC", "BLGSP_Bcell_UNC", "BL_Pediatric") | (sample_id == "06-29223T"))

    }else if(case_set == "BL-DLBCL-manuscript"){
      adult_bl_manuscript_samples = data.table::fread("/projects/rmorin/projects/gambl-repos/gambl-kdreval/data/metadata/BLGSP--DLBCL-case-set.tsv") %>%
        pull(Tumor_Sample_Barcode)

      all_meta = all_meta %>%
      dplyr::filter(sample_id %in% adult_bl_manuscript_samples)

    }else if(case_set == "FL-DLBCL-all"){
      fl_dlbcl_all_samples = data.table::fread("/projects/rmorin/projects/gambl-repos/gambl-kdreval/data/metadata/FL--DLBCL--all-case-set.tsv") %>%
        pull(Tumor_Sample_Barcode)

      all_meta = all_meta %>%
      dplyr::filter(sample_id %in% fl_dlbcl_all_samples)

    }else if(case_set == "GAMBL-all"){

      #get all GAMBL but remove FFPE benchmarking cases and ctDNA
      all_meta = all_meta %>%
      dplyr::filter(!cohort %in% c("FFPE_Benchmarking", "DLBCL_ctDNA"))
    }else if(case_set %in% colnames(full_case_set)) {
      # ensure consistent column naming
      full_case_set =
        full_case_set %>% rename_at(vars(matches(
          "sample_id", ignore.case = TRUE
        )),
        ~ "Tumor_Sample_Barcode")

      # get case set as defined in the file
      this_subset_samples =
        full_case_set %>%
        dplyr::filter(!!sym(case_set) == 1) %>%
        pull(Tumor_Sample_Barcode)

      all_meta = all_meta %>%
        dplyr::filter(sample_id %in% this_subset_samples)

    }else{
      message(paste("case set", case_set, "not available"))
      return()
    }
  }

  #add some derivative columns that simplify and consolidate some of the others (DLBCL-specific)
  #all_meta = all_meta %>% dplyr::mutate(lymphgen = case_when(
  #  pathology != "DLBCL" ~ pathology,
  #  str_detect(lymphgen_cnv_noA53,"/") ~ "COMPOSITE",
  #  TRUE ~ lymphgen_cnv_noA53
  #))

  all_meta = GAMBLR::tidy_lymphgen(all_meta,
              lymphgen_column_in = "lymphgen_cnv_noA53",
              lymphgen_column_out = "lymphgen",
              relevel=TRUE)

  all_meta = GAMBLR::collate_lymphgen(all_meta, verbose=FALSE)

  # "catchall" pathology for those that need review
  all_meta = all_meta %>%
    mutate(pathology = ifelse(nchar(pathology) > 15, "OTHER", pathology))

  all_meta = mutate(all_meta, Tumor_Sample_Barcode = sample_id) #duplicate for convenience
  all_meta = all_meta %>%
    dplyr::mutate(consensus_coo_dhitsig = case_when(pathology != "DLBCL" ~ pathology,
                                                    COO_consensus == "ABC" ~ COO_consensus,
                                                    DLBCL90_dhitsig_call == "POS" ~ "DHITsigPos",
                                                    DLBCL90_dhitsig_call == "NEG" ~ "DHITsigNeg",
                                                    DHITsig_PRPS_class == "DHITsigPos" ~ "DHITsigPos",
                                                    DHITsig_PRPS_class == "DHITsig+" ~ "DHITsigPos",
                                                    DHITsig_PRPS_class == "DHITsigNeg" ~ "DHITsigNeg",
                                                    DHITsig_PRPS_class == "DHITsig-" ~ "DHITsigNeg",
                                                    DHITsig_PRPS_class == "UNCLASS" ~ "DHITsigPos",
                                                    TRUE ~ "NA"))

  all_meta = all_meta %>%
    dplyr::mutate(DHITsig_consensus = case_when(DHITsig_consensus == "NEG" ~ "DHITsigNeg",
                                                DHITsig_consensus == "POS" ~ "DHITsigPos",
                                                DHITsig_consensus == "UNCLASS" ~ "DHITsig-IND",
                                                DHITsig_consensus == "DHITsigNeg" ~ DHITsig_consensus,
                                                DHITsig_consensus == "DHITsigPos" ~ DHITsig_consensus,
                                                TRUE ~ "NA"))

  #assign a rank to each pathology for consistent and sensible ordering
  all_meta = all_meta %>%
    dplyr::mutate(pathology_rank = case_when(pathology == "B-ALL" ~ 0,
                                             pathology == "SCBC" ~ 2,
                                             pathology == "CLL" ~ 3,
                                             pathology == "MCL" ~ 4,
                                             pathology == "BL" ~ 7,
                                             pathology == "DLBCL-BL-like" ~ 10,
                                             pathology == "HGBL" ~ 11,
                                             pathology == "COMFL" ~ 13,
                                             pathology == "FL" ~ 15,
                                             pathology == "DLBCL" ~ 19,
                                             pathology == "PBL" ~ 27,
                                             pathology == "B-cell unclassified" ~ 29,
                                             pathology == "MM" ~ 33,
                                             TRUE ~ 35))

  all_meta = all_meta %>%
    dplyr::mutate(lymphgen_rank = case_when(pathology != "DLBCL" ~ pathology_rank,
                                            lymphgen == "Other" ~ 16,
                                            lymphgen == "COMPOSITE" ~ 17,
                                            lymphgen == "N1" ~ 18,
                                            lymphgen == "EZB" ~ 19,
                                            lymphgen == "ST2" ~ 20,
                                            lymphgen == "BN2" ~ 21,
                                            lymphgen == "MCD" ~ 22,
                                            TRUE ~ 50))
  if(with_outcomes){

    all_meta = left_join(all_meta, outcome_table, by = "patient_id") %>%
      mutate(age_group = case_when(cohort == "BL_Adult"~"Adult_BL", cohort == "BL_Pediatric" | cohort == "BL_ICGC" ~ "BL_Pediatric", TRUE ~ "Other"))

  }
  # take one row per sample_id using seq_type_priority
  if(seq_type_priority=="genome"){
    all_meta = all_meta %>%
      arrange(sample_id,seq_type) %>%
      group_by(sample_id) %>%
      slice_tail() %>%
      ungroup()
  }
  if(seq_type_priority=="capture"){
    all_meta = all_meta %>%
      arrange(sample_id,seq_type) %>%
      group_by(sample_id) %>%
      slice_head() %>%
      ungroup()
  }
  return(all_meta)
}

add_prps_result = function(incoming_metadata){
  prps_res = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/derived_and_curated_metadata/outputs/BL_dhitsig_PRPS.tsv"))
  colnames(prps_res)[1] = "sample_id"
  prps_res = dplyr::select(prps_res, sample_id, PRPS_score, PRPS_class)

  #need to associate each sample with a patient ID then annotate the metadata based on patient ID
  patient_meta_g = get_gambl_metadata(seq_type_filter = "genome") %>%
    dplyr::select(sample_id, patient_id)

  patient_meta_r = get_gambl_metadata(seq_type_filter = "mrna") %>%
    dplyr::select(sample_id, patient_id)

  patient_meta = bind_rows(patient_meta_g, patient_meta_r)
}


#' INTERNAL FUNCTION called by get_gambl_metadata, not meant for out-of-package usage.
#' Layer on ICGC metadata from a supplemental table to fill in missing COO.
#'
#' @param incoming_metadata A metadata table (probably output from get_gambl_metadata).
#'
#' @return Meta data with layered information (ICGC).
#'
#' @examples
#' icgc_metadata = add_icgc_metadata(incoming_metadata = my_meta)
#'
add_icgc_metadata = function(incoming_metadata){
  repo_base = config::get("repo_base")
  icgc_publ_file = paste0(repo_base,"data/metadata/raw_metadata/MALY_DE_tableS1.csv")
  icgc_publ = suppressMessages(suppressWarnings(read_csv(icgc_publ_file)))
  icgc_publ = icgc_publ[,c(1:20)]
  #fix commas as decimals
  icgc_publ = mutate(icgc_publ, purity = str_replace(purity, ",", "."))
  icgc_publ = mutate(icgc_publ, sex = str_to_upper(sex))

  icgc_raw_path = paste0(repo_base,"data/metadata/raw_metadata/ICGC_MALY_seq_md.tsv")
  icgc_raw = suppressMessages(read_tsv(icgc_raw_path))

  icgc_raw = icgc_raw %>%
    dplyr::select(-compression, -bam_available, -read_length, -time_point, -unix_group, -ffpe_or_frozen, -link_name)  %>%
    dplyr::filter(tissue_status %in% c("tumor", "tumour"))

  icgc_all = left_join(icgc_raw, icgc_publ,by = "ICGC_ID") %>%
    dplyr::select(-tissue_status, -seq_type, -protocol, -seq_source_type, -data_path, -genome_build, -RNA_available) %>%
    dplyr::select(sample_id, ICGC_ID, pathology.x, pathology.y, COO, molecular_BL, MYC_sv, BCL2_sv, BCL6_sv) %>%
    dplyr::rename("ICGC_MYC_sv" = "MYC_sv") %>%
    dplyr::rename("ICGC_BCL2_sv" = "BCL2_sv") %>%
    dplyr::rename("ICGC_BCL6_sv" = "BCL6_sv") %>%
    dplyr::rename("detailed_pathology" = "pathology.x") %>%
    dplyr::rename("ICGC_PATH" = "pathology.y")

  #join with all metadata to fill in blanks
  #all_meta=get_gambl_metadata()
  rejoined = left_join(incoming_metadata, icgc_all,by = "sample_id") %>%
    mutate(COO_final = case_when(!is.na(COO_consensus) ~ COO_consensus, COO != "n.a." & COO != "TypeIII" ~ COO, TRUE ~ "NA")) %>%
    dplyr::select(-COO)
  return(rejoined)
}


#' INTERNAL FUNCTION called by get_gambl_metadata, not meant for out-of-package usage.
#' Get the patient-centric clinical metadata.
#'
#' @param patient_ids Vector of patient IDs.
#' @param time_unit Return follow-up times in one of three time units: year, month or day.
#' @param censor_cbioportal Optionally request the censoring to be encoded in the specific style required by cBioPortal.
#' @param complete_missing Optionally fill in any gaps to ensure we have values for every patient (censor at 0 if missing).
#' @param from_flatfile Optionally set to FALSE to use the database to get the survival data.
#'
#' @return Data frame with one row for each patient_id.
#' @import tidyverse RMariaDB DBI dbplyr
#'
#' @examples
#' outcome_df = get_gambl_outcomes()
#'
get_gambl_outcomes = function(patient_ids,
                              time_unit = "year",
                              censor_cbioportal = FALSE,
                              complete_missing = FALSE,
                              from_flatfile = TRUE){

  if(from_flatfile){
    outcome_flatfile = paste0(config::get("repo_base"), config::get("table_flatfiles")$outcomes)
    all_outcome = suppressMessages(read_tsv(outcome_flatfile))

  }else{
    db = config::get("database_name")
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
    all_outcome = dplyr::tbl(con, "outcome_metadata") %>%
      as.data.frame()

    DBI::dbDisconnect(con)
  }
  if(!missing(patient_ids)){
    all_outcome = all_outcome %>%
      dplyr::filter(patient_id %in% patient_ids)

    if(complete_missing){
      #add NA values and censored outcomes for all missing patient_ids
      all_outcome = all_outcome %>%
        complete(patient_id = patient_ids, fill = list(OS_YEARS = 0, PFS_years = 0, TTP_YEARS = 0, DSS_YEARS = 0, CODE_OS = 0, CODE_PFS = 0, CODE_DSS = 0, CODE_TTP = 0))
    }
  }
  if(time_unit == "month"){
    all_outcome = all_outcome %>%
      mutate(OS_MONTHS = OS_YEARS * 12)

    all_outcome = all_outcome %>%
      mutate(PFS_MONTHS = PFS_YEARS * 12)

    all_outcome = all_outcome %>%
      mutate(TTP_MONTHS = TTP_YEARS * 12)

    all_outcome = all_outcome %>%
      mutate(DSS_MONTHS = DSS_YEARS * 12)

    all_outcome = all_outcome %>%
      dplyr::select(-c("OS_YEARS", "PFS_YEARS", "TTP_YEARS", "DSS_YEARS"))

  }else if(time_unit == "day"){
    all_outcome = all_outcome %>%
      mutate(OS_DAYS = OS_YEARS * 365)

    all_outcome = all_outcome %>%
      mutate(PFS_DAYS = PFS_YEARS * 365)

    all_outcome = all_outcome %>%
      mutate(TTP_DAYS = TTP_YEARS * 365)

    all_outcome = all_outcome %>%
      mutate(DSS_DAYS = DSS_YEARS * 365)

    all_outcome = all_outcome %>%
      dplyr::select(-c("OS_YEARS", "PFS_YEARS", "TTP_YEARS", "DSS_YEARS"))
  }

  #if necessary, convert the censoring into the cBioPortal format for OS and PFS
  if(censor_cbioportal){
    all_outcome$OS_STATUS = as.character(all_outcome$CODE_OS)
    all_outcome = all_outcome %>%
      mutate(OS_STATUS = case_when(OS_STATUS == "0" ~ "0:LIVING", OS_STATUS == "1"~"1:DECEASED"))

    all_outcome$DFS_STATUS = as.character(all_outcome$CODE_PFS)
    all_outcome = all_outcome %>%
      mutate(DFS_STATUS = case_when(DFS_STATUS == "0" ~ "0:DiseaseFree", DFS_STATUS == "1"~"1:Recurred/Progressed"))

    all_outcome = all_outcome %>%
      mutate(all_outcome, DFS_MONTHS = PFS_MONTHS)
  }
  all_outcome = all_outcome %>%
    mutate(is_adult = ifelse(age < 20, "Pediatric", "Adult"))

  return(all_outcome)
}


#' Retrieve Combined Manta and GRIDSS-derived SVs from a flatfile and filter.
#'
#' The bedpe files used as input to this function were pre-filtered for a minimum VAF of 0.05, and SVs affecting.
#' common translocation regions (BCL2, BCL6, MYC, CCND1) were whitelisted (e.g. no VAF filter applied).
#' Therefore if you wish to post-filter the SVs we recommend doing so carefully after loading this data frame.
#' Further, the input bedpe file is annotated with oncogenes and superenhancers from naive and germinal centre B-cells.
#' You can subset to events affecting certain loci using the "oncogenes" argument.
#'
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Recommended: 0. (default: 0)
#' @param sample_ids A character vector of tumour sample IDs you wish to retrieve SVs for.
#' @param with_chr_prefix Prepend all chromosome names with chr (required by some downstream analyses)
#' @param projection The projection genome build. Currently only grch37 is supported but hg38 should be easy to add.
#' @param oncogenes A character vector of genes commonly involved in translocations. Possible values: CCND1, CIITA, SOCS1, BCL2, RFTN1, BCL6, MYC, PAX5.
#'
#'
#' @return A data frame in a bedpe-like format with additional columns that allow filtering of high-confidence SVs.
#' @export
#' @import DBI RMariaDB tidyverse dbplyr
#'
#' @examples
#' get_combined_sv(oncogenes = c("MYC", "BCL2", "BCL6"))
#'
get_combined_sv = function(min_vaf = 0,
                           sample_ids,
                           with_chr_prefix = FALSE,
                           projection = "grch37",
                           oncogenes){

  if(projection != "grch37"){
    message("Currently, only grch37 is supported")

  return()
  }

  base_path = config::get("project_base")
  sv_file = config::get()$results_flatfiles$sv_combined$icgc_dart
  if(projection == "hg38"){
    sv_file = str_replace(sv_file, "--grch37", "--hg38")
  }
  sv_file = paste0(base_path, sv_file)
  permissions = file.access(sv_file, 4)
  if(permissions == - 1){
    sv_file = config::get()$results_flatfiles$sv_combined$gambl
    sv_file = paste0(base_path, sv_file)
  }
  all_sv = read_tsv(sv_file, col_types = "cnncnncnccccnnccncn") %>%
    dplyr::rename(c("VAF_tumour" = "VAF")) %>%
    dplyr::filter(VAF_tumour >= min_vaf)

  if(!missing(sample_ids)){
    all_sv = all_sv %>%
      dplyr::filter(tumour_sample_id %in% sample_ids)
  }

  if(!missing(oncogenes)){
    all_sv = all_sv %>%
      dplyr::filter(ANNOTATION_A %in% oncogenes | ANNOTATION_B %in% oncogenes)
  }

  if(with_chr_prefix){
    #add chr prefix only if it's missing
    all_sv = all_sv %>%
      dplyr::mutate(CHROM_A = case_when(str_detect(CHROM_A, "chr") ~ CHROM_A, TRUE ~ paste0("chr", CHROM_A)))

    all_sv = all_sv %>%
      dplyr::mutate(CHROM_B = case_when(str_detect(CHROM_B, "chr") ~ CHROM_B, TRUE ~ paste0("chr", CHROM_B)))
  }
  return(all_sv)
}


#' Retrieve Manta SVs from the database and filter.
#'
#' @param min_vaf The minimum tumour VAF for a SV to be returned.
#' @param min_score The lowest Manta somatic score for a SV to be returned.
#' @param pass If set to TRUE, include SVs that are annotated with PASS in FILTER column. Default is TRUE.
#' @param pair_status Use to restrict results (if desired) to matched or unmatched results (default is to return all).
#' @param sample_id Filter on specific sample IDs in tumour_sample_id column.
#' @param chromosome The chromosome you are restricting to.
#' @param qstart Query start coordinate of the range you are restricting to.
#' @param qend Query end coordinate of the range you are restricting to.
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param with_chr_prefix Prepend all chromosome names with chr (required by some downstream analyses).
#' @param from_flatfile Set to FALSE by default.
#' @param projection The projection genome build. Currently only grch37 is supported.
#'
#' @return A data frame in a bedpe-like format with additional columns that allow filtering of high-confidence SVs.
#' @export
#' @import DBI RMariaDB tidyverse dbplyr
#'
#' @examples
#' #lazily get every SV in the table with default quality filters
#' all_sv = get_manta_sv()
#' #get all SVs for a single sample
#' some_sv = get_manta_sv(sample_id="94-15772_tumorA")
#' #get the SVs in a region around MYC
#' myc_locus_sv = get_manta_sv(region="8:128723128-128774067").
#'
get_manta_sv = function(min_vaf = 0.1,
                        min_score = 40,
                        pass = TRUE,
                        pairing_status,
                        sample_id,
                        chromosome,
                        qstart,
                        qend,
                        region,
                        with_chr_prefix = FALSE,
                        from_flatfile = TRUE,
                        projection = "grch37"){

  if(projection != "grch37"){
    message("Currently, only grch37 is supported")
    return()
  }

  db = config::get("database_name")
  table_name = config::get("results_tables")$sv
    if(!missing(region)){
      region = gsub(",", "", region)
      #format is chr6:37060224-37151701
      split_chunks = unlist(strsplit(region, ":"))
      chromosome = split_chunks[1]
      startend = unlist(strsplit(split_chunks[2], "-"))
      qstart = startend[1]
      qend = startend[2]
    }

  #this table stores chromosomes with un-prefixed names. Convert to prefixed chromosome if necessary
  if(from_flatfile){
    sv_file = get_merged_result(tool_name = "manta", projection = projection)
    all_sv = read_tsv(sv_file, col_types = "cnncnncnccccnnnnccc", col_names = cnames)
  }else{
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
    all_sv = dplyr::tbl(con, table_name) %>%
      as.data.frame()
  }
  if(!missing(region) || !missing(chromosome)){
    suppressWarnings({
      if(grepl("chr",chromosome)){
        chromosome = gsub("chr", "", chromosome)
      }
    })
    all_sv = all_sv %>%
      dplyr::filter((CHROM_A == chromosome & START_A >= qstart & START_A <= qend) | (CHROM_B == chromosome & START_B >= qstart & START_B <= qend)) %>%
      dplyr::filter(VAF_tumour >= min_vaf & SOMATIC_SCORE >= min_score)
  }else{
    all_sv = all_sv %>%
      dplyr::filter(VAF_tumour >= min_vaf & SOMATIC_SCORE >= min_score)
  }
  if(pass){
    all_sv = all_sv %>%
      dplyr::filter(FILTER == "PASS")
  }
  if(!missing(pairing_status)){
    all_sv = all_sv %>%
      dplyr::filter(pair_status == pairing_status)
  }
  if(!missing(sample_id)){
    all_sv = all_sv %>%
      dplyr::filter(tumour_sample_id == sample_id)
  }
  all_sv = as.data.frame(all_sv)
  if(with_chr_prefix){
    #add chr prefix only if it's missing
    all_sv = all_sv %>%
      dplyr::mutate(CHROM_A = case_when(str_detect(CHROM_A, "chr") ~ CHROM_A, TRUE ~ paste0("chr", CHROM_A)))

    all_sv = all_sv %>%
      dplyr::mutate(CHROM_B = case_when(str_detect(CHROM_B, "chr") ~ CHROM_B, TRUE ~ paste0("chr", CHROM_B)))
  }
  if(!from_flatfile){
    DBI::dbDisconnect(con)
  }
  return(all_sv)
}


#' Get a copy number matrix for all samples based on segmented data in database.
#'
#' @param regions_list A list of regions in the format chrom:start-end.
#' @param regions_bed A bed file with one row for each region you want to determine the CN state from.
#' @param region_names Subset CN states on specific regions (gene symbols e.g FCGR2B).
#' @param all_cytobands Include all cytobands, default is set to FALSE. Currently only supports hg19.
#' @param use_cytoband_name Use cytoband names instad of region name, e.g p36.33.
#'
#' @return Copy number matrix.
#' @import tidyverse circlize
#' @export
#'
#' @examples
#' #basic usage, generic lymphoma gene list
#' cn_matrix = get_cn_states(regions_bed=grch37_lymphoma_genes_bed)
#' single_gene_cn = get_cn_states(regions_list=c(this_region), region_names = c("FCGR2B"))
#'
get_cn_states = function(regions_list,
                         regions_bed,
                         region_names,
                         all_cytobands = FALSE,
                         use_cytoband_name = FALSE){
  if(all_cytobands){
    message("Currently, only grch37 is supported")
  }
  #retrieve the CN value for this region for every segment that overlaps it
  bed2region=function(x){
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }
  if(all_cytobands){
    message("Cytobands are in respect to hg19. This will take awhile but it does work, trust me!")
    regions_bed = circlize::read.cytoband(species = "hg19")$df
    colnames(regions_bed) = c("chromosome_name", "start_position", "end_position", "name", "dunno")
    if(use_cytoband_name){
      regions_bed = mutate(regions_bed, region_name = paste0(str_remove(chromosome_name, pattern = "chr"), name))
      region_names = pull(regions_bed, region_name)
    }else{
      region_names = pull(regions_bed, region_name)
    }
    regions = apply(regions_bed, 1, bed2region)
    #use the cytobands from the circlize package (currently hg19 but can extend to hg38 once GAMBLR handles it) Has this been updated?
  }else if(missing(regions_list)){
    if(!missing(regions_bed)){
      regions = apply(regions_bed, 1, bed2region)
    }else{
      warning("You must supply either regions_list or regions_df")
    }
  }else{
    regions = regions_list
  }
  region_segs = lapply(regions,function(x){get_cn_segments(region = x, streamlined = TRUE)})
  if(missing(region_names) & !use_cytoband_name){
    region_names = regions
  }
  tibbled_data = tibble(region_segs, region_name = region_names)
  unnested_df = tibbled_data %>%
    unnest_longer(region_segs)

  seg_df = data.frame(ID = unnested_df$region_segs$ID, CN = unnested_df$region_segs$CN,region_name = unnested_df$region_name)
  #arbitrarily take the first segment for each region/ID combination
  seg_df = seg_df %>%
    dplyr::group_by(ID, region_name) %>%
    dplyr::slice(1) %>%
    dplyr::rename("sample_id" = "ID")

  #fill in any sample/region combinations with missing data as diploid
  meta_arranged = get_gambl_metadata() %>%
    dplyr::select(sample_id, pathology, lymphgen) %>%
    arrange(pathology, lymphgen)

  eg = expand_grid(sample_id = pull(meta_arranged, sample_id), region_name = as.character(unique(seg_df$region_name)))
  all_cn = left_join(eg, seg_df, by = c("sample_id" = "sample_id", "region_name" = "region_name")) %>%
    mutate(CN = replace_na(CN, 2))

  cn_matrix = pivot_wider(all_cn, id_cols = "sample_id", names_from = "region_name", values_from = "CN") %>%
    column_to_rownames("sample_id")

  names(cn_matrix) = region_names

  #order the regions the same way the user provided them for convenience
  return(cn_matrix)
}


#' Get all segments for a single (or multiple) sample_id(s).
#'
#' @param this_sample_id Optional argument, single sample_id for the sample to retrieve segments for.
#' @param multiple_samples Set to TRUE to return cn segments for multiple samples (list) pf samples to be specified in samples_list parameter.
#' @param samples_list Optional vector of type character with one sample per row, rewuired if multiple_samples is set to TRUE.
#' @param with_chr_prefix Set to TRUE to add a chr prefix to chromosome names.
#' @param streamlined Return a minimal output rather than full details.
#'
#' @return A list of segments for a specific or multiple sample ID(s)
#' @export
#'
#' @examples
#' Return cn segments for one sample:
#' sample_cn_seg = get_sample_cn_segments(this_sample_id = "some-sample-id", multiple_samples = FALSE)
#'
#' Return cn segments for multiple samples (provided as list):
#' samples = get_sample_cn_segments(multiple_samples = TRUE, sample_list = c("some_sample", "another_sample"))
#'
#' Return cn segments for multiple samples (read csv with one sample per line):
#' sample_list = readLines("../samples-test.csv")
#' multiple_samples = get_sample_cn_segments(multiple_samples = TRUE, sample_list = sample_list)
#'
get_sample_cn_segments = function(this_sample_id,
                                  multiple_samples = FALSE,
                                  sample_list,
                                  from_flatfile = TRUE,
                                  projection = "grch37",
                                  with_chr_prefix = FALSE,
                                  streamlined = FALSE){
  if(from_flatfile){
    seq_type = "genome"
    cnv_flatfile_template = config::get("results_flatfiles")$cnv_combined$icgc_dart
    cnv_path =  glue::glue(cnv_flatfile_template)
    full_cnv_path =  paste0(config::get("project_base",config="default"), cnv_path)
    local_full_cnv_path =  paste0(config::get("project_base"), cnv_path)
    if(file.exists(local_full_cnv_path)){
      full_cnv_path = local_full_cnv_path
      #use local file when available
    }
    # check permissions to ICGC data
    permissions = file.access(full_cnv_path, 4)
    if (permissions == -1) {
      message("restricting to non-ICGC data")
      cnv_flatfile_template = config::get("results_flatfiles")$cnv_combined$gambl
      cnv_path =  glue::glue(cnv_flatfile_template)
      full_cnv_path =  paste0(config::get("project_base"), cnv_path)
    }

    all_segs = read_tsv(full_cnv_path)
    if (!missing(this_sample_id) & !multiple_samples) {
      all_segs = dplyr::filter(all_segs, ID %in% this_sample_id)
    } else if (!missing(sample_list)) {
      all_segs = dplyr::filter(all_segs, ID %in% sample_list)
    }

  } else {

    if(!missing(this_sample_id) & !multiple_samples){
      sample_status = get_gambl_metadata() %>%
        dplyr::filter(sample_id == this_sample_id) %>%
        pull(pairing_status)

      db = config::get("database_name")
      table_name = config::get("results_tables")$copy_number
      table_name_unmatched = config::get("results_tables")$copy_number_unmatched
      con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

      all_segs_matched = dplyr::tbl(con, table_name) %>%
        dplyr::filter(ID == this_sample_id) %>%
        as.data.frame() %>%
        dplyr::mutate(method = "battenberg")

      all_segs_unmatched = dplyr::tbl(con, table_name_unmatched) %>%
        dplyr::filter(ID == this_sample_id) %>%
        as.data.frame() %>%
        dplyr::filter(! ID %in% all_segs_matched$ID) %>%
        dplyr::mutate(method = "controlfreec")
    }

    if(multiple_samples & missing(this_sample_id)){
      sample_status = get_gambl_metadata() %>%
        dplyr::filter(sample_id %in% sample_list) %>%
        pull(pairing_status)

      db = config::get("database_name")
      table_name = config::get("results_tables")$copy_number
      table_name_unmatched = config::get("results_tables")$copy_number_unmatched
      con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

      all_segs_matched = dplyr::tbl(con, table_name) %>%
        dplyr::filter(ID %in% sample_list) %>%
        as.data.frame() %>%
        dplyr::mutate(method = "battenberg")

      all_segs_unmatched = dplyr::tbl(con, table_name_unmatched) %>%
        dplyr::filter(ID %in% sample_list) %>%
        as.data.frame() %>%
        dplyr::filter(! ID %in% all_segs_matched$ID) %>%
        dplyr::mutate(method = "controlfreec")
    }

    all_segs = rbind(all_segs_matched, all_segs_unmatched)

  }


  all_segs = dplyr::mutate(all_segs, CN = round(2*2^log.ratio))

  if(!with_chr_prefix){all_segs = all_segs %>% dplyr::mutate(chrom = gsub("chr", "", chrom))}
  if(streamlined){all_segs = dplyr::select(all_segs, ID, CN)}

  return(all_segs)
}


#' Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range.
#'
#' @param chromosome The chromosome you are restricting to.
#' @param qstart Start coordinate of the range you are restricting to.
#' @param qend End coordinate of the range you are restricting to.
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param with_chr_prefix Prepend all chromosome names with chr (required by some downstream analyses).
#' @param streamlined Return a basic rather than full MAF format.
#' @param from_flatfile Set to FALSE by default.
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
#' #basic usage
#' my_segments = get_cn_segments(region="chr8:128,723,128-128,774,067")
#' #specifying chromosome, start and end individually
#' my_segments = get_cn_segments(chromosome="8",qstart=128723128,qend=128774067)
#' #Asking for chromosome names to have a chr prefix (default is un-prefixed)
#' prefixed_segments = get_cn_segments(chromosome ="12",qstart = 122456912, qend = 122464036, with_chr_prefix = TRUE)
#'
get_cn_segments = function(chromosome = "",
                           qstart,
                           qend,
                           region,
                           with_chr_prefix = FALSE,
                           streamlined = FALSE,
                           from_flatfile = FALSE,
                           ssh_session){

  db = config::get("database_name")
  table_name = config::get("results_tables")$copy_number
  table_name_unmatched = config::get("results_tables")$copy_number_unmatched
  if(!missing(region)){
    region = gsub(",", "", region)
    #format is chr6:37060224-37151701
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }
  if(grepl("chr", chromosome)){
  }else{
    chromosome = paste0("chr", chromosome)
  }

  #chr prefix the query chromosome to match how it's stored in the table.
  #This isn't yet standardized in the db so it's just a workaround "for now".
  if(from_flatfile){
    base_dir = config::get(config="default")$project_base

    unmatched_path = config::get()$results_directories$controlfreec

    #separated by which genome the sample was aligned to
    unmatched_hg38_path = paste0(unmatched_path, "from--genome--hg38/")
    unmatched_grch37_path = paste0(unmatched_path, "from--genome--grch37/")
  }else{
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

    #remove the prefix if this is false (or leave as is otherwise)
    if(missing(qstart) & missing(region)){
      all_segs_matched = dplyr::tbl(con, table_name) %>%
        as.data.frame()

      all_segs_unmatched = dplyr::tbl(con, table_name_unmatched) %>%
        as.data.frame()
    }else{

      #TODO improve this query to allow for partial overlaps, create Issue on Github?
      all_segs_matched = dplyr::tbl(con, table_name) %>%
        dplyr::filter((chrom == chromosome & start <= qstart & end >= qend) | (chrom == chromosome & start >= qstart & end <= qend)) %>%
        as.data.frame() %>%
        dplyr::mutate(method = "battenberg")

      # get controlfreec segments for samples with missing battenberg results like unpaired
      all_segs_unmatched = dplyr::tbl(con, table_name_unmatched) %>%
        dplyr::filter((chrom == chromosome & start <= qstart & end >= qend) | (chrom == chromosome & start >= qstart & end <= qend)) %>%
        as.data.frame() %>%
        dplyr::filter(! ID %in% all_segs_matched$ID)  %>%
        dplyr::mutate(method = "controlfreec")

      DBI::dbDisconnect(con)
    }
  }
  all_segs = rbind(all_segs_matched, all_segs_unmatched)
  all_segs = dplyr::mutate(all_segs, CN = round(2*2^log.ratio))
  if(! with_chr_prefix){
    all_segs = all_segs %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom))
  }
  if(streamlined){
    all_segs = dplyr::select(all_segs, ID, CN)
  }
  return(all_segs)
}


#' INTERNAL FUNCTION, not meant for out-of-package usage.
#' Housekeeping function to add results to a table.
#'
#' @param table_name The name of the database table to update/populate.
#' @param data_df A dataframe of values to load into the table.
#'
#' @return Table
#'
#' @examples
#' table_up = append_to_table("my_table", "my_df")
#'
append_to_table = function(table_name,
                           data_df){

  db = config::get("database_name")
  con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  dbWriteTable(con, table_name, table_data, append = TRUE)
}


#' Prepare a matrix with one row per sample and one column per region using a set of hypermutated regions.
#' Values are the number of mutations in that patient in the region.
#'
#' @param regions_bed A bed file with one row for each region.
#' @param maf_data Optionally provide a data frame in the MAF format, otherwise the database will be used.
#' @param sample_metadata This is used to complete your matrix. All GAMBL samples will be used by default. Provide a data frame with at least sample_id for all samples if you are using non-GAMBL data.
#' @param use_name_column Set this to true to force the function to use the value in column "name" to name each feature in the output.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details).
#' @param allow_clustered Set to TRUE to utilize the latest SLMS-3 variant calls that allow clustered variants.
#'
#' @return Matrix
#' @export
#'
#' @examples
#' matrix = get_ashm_count_matrix(regions_bed = "my_bed.bed", sample_metadata = "GAMBL-metadata")
#'
get_ashm_count_matrix = function(regions_bed,
                                 maf_data,
                                 sample_metadata,
                                 use_name_column = FALSE,
                                 from_indexed_flatfile = TRUE){

  if(missing(regions_bed)){
    regions_bed = grch37_ashm_regions
  }
  ashm_maf = suppressMessages(get_ssm_by_regions(regions_bed = regions_bed,
                                basic_columns = TRUE,
                                maf_data = maf_data,
                                use_name_column = use_name_column,
                                from_indexed_flatfile = from_indexed_flatfile))

  ashm_counted = ashm_maf %>%
    group_by(Tumor_Sample_Barcode, region_name) %>%
    tally()

  colnames(ashm_counted)[1] = "sample_id"

  if(missing(sample_metadata)){
    all_meta = get_gambl_metadata() %>%
      dplyr::select(sample_id)
  }else{
    all_meta = dplyr::select(sample_metadata, sample_id)
  }

  #fill out all combinations so we can get the cases with zero mutations
  eg = expand_grid(sample_id = pull(all_meta, sample_id), region_name = unique(ashm_counted$region_name))
  all_counts = left_join(eg, ashm_counted) %>%
    mutate(n = replace_na(n, 0)) %>%
    unique() #not sure where the duplicates are coming from but its annoying

  all_counts_wide = pivot_wider(all_counts, id_cols = sample_id, names_from = region_name, values_from = n) %>%
    column_to_rownames(var = "sample_id")

  return(all_counts_wide)
}


#' Efficiently retrieve all mutations across a range of genomic regions.
#'
#' @param regions_list Either provide a vector of regions in the chr:start-end format OR.
#' @param regions_bed Better yet, provide a bed file with the coordinates you want to retrieve.
#' @param streamlined Set to TRUE to return only the following columns; start, sample_id and region_name.
#' @param basic_columns Set to FALSE to override the default behavior of returning only the following columns; chromosome, start, end and sample_id.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a list of indexes (integer) or a list of characters.
#' @param return_cols If set to TRUE, a vector with all available column names will be returned. Default is FALSE.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details).
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF
#' @param seq_type The seq_type you want back, default is genome.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param allow_clustered Logical parameter indicating whether to use SLMS-3 results with clustered events. Default is FALSE.
#'
#' @return Returns a data frame of variants in MAF-like format.
#' @export
#'
#' @examples
#' #basic usage, adding custom names from bundled ashm data frame
#' regions_bed = grch37_ashm_regions %>% mutate(name = paste(gene, region, sep = "_"))
#' ashm_maf = get_ssm_by_regions(regions_bed = regions_bed)
#'
get_ssm_by_regions = function(regions_list,
                              regions_bed,
                              streamlined = FALSE,
                              basic_columns = TRUE,
                              maf_cols = NULL,
                              return_cols = FALSE,
                              use_name_column = FALSE,
                              maf_data = maf_data,
                              from_indexed_flatfile = TRUE,
                              mode = "slms-3",
                              augmented = TRUE,
                              seq_type = "genome",
                              projection = "grch37",
                              min_read_support = 4){

  bed2region = function(x){
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }

  if(missing(regions_list)){
    if(!missing(regions_bed)){
      regions = apply(regions_bed, 1, bed2region)
    }else{
      warning("You must supply either regions_list or regions_df")
    }
  }
  if(missing(maf_data)){
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x,
                                                                streamlined = streamlined,
                                                                basic_columns = TRUE,
                                                                from_indexed_flatfile = from_indexed_flatfile,
                                                                mode = mode,
                                                                augmented = augmented,
                                                                seq_type = seq_type,
                                                                projection = projection)})
  }else{
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x,
                                                                streamlined = streamlined,
                                                                maf_data = maf_data,
                                                                from_indexed_flatfile = from_indexed_flatfile,
                                                                mode = mode)})
  }

  if(!use_name_column){
    rn = regions
  }else{
    rn = regions_bed[["name"]]
  }

  #get region names
  tibbled_data = tibble(region_mafs, region_name = rn)
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)

  #subset on region names and rename column before cbind with full maf
  r_names = as.data.frame(unnested_df$region_name)
  colnames(r_names)[1] <- "region_name"

  #unnest maf
  unnested_maf = bind_rows(region_mafs)

  #join region names df with unnested maf
  maf = cbind(r_names, unnested_maf)

  if(streamlined){
    maf = maf %>%
      dplyr::select(Start_Position, Tumor_Sample_Barcode)
  }else if(basic_columns){
    maf = maf[, c(1:46)]
  }

  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    maf = dplyr::select(maf, all_of(maf_cols))
  }

  #print all available columns
  if(!basic_columns && return_cols){
    print(colnames(maf))
  }

  return(maf)
}


#' Retrieve all SSMs from the GAMBL database within a single genomic coordinate range.
#'
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param qstart Query start coordinate of the range you are restricting to.
#' @param qend Query end coordinate of the range you are restricting to.
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param basic_columns Set to FALSE to override the default behavior of returning only the first 45 columns of MAF data.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a list of indexes (integer) or a list of characters.
#' @param return_cols If set to TRUE, a vector with all available column names will be returned. Default is FALSE.
#' @param streamlined Return a basic rather than full MAF format, default is FALSE.
#' @param maf_data Parameter description.
#' @param seq_type The seq_type you want back, default is genome.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details).
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF .
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' #basic usage
#' my_mutations = get_ssm_by_region(region = "chr8:128,723,128-128,774,067")
#' #specifying chromosome, start and end individually
#' my_mutations = get_ssm_by_region(chromosome = "8", qstart = 128723128, qend = 128774067)
#'
get_ssm_by_region = function(chromosome,
                             qstart,
                             qend,
                             region = "",
                             basic_columns = TRUE,
                             maf_cols = NULL,
                             return_cols = FALSE,
                             streamlined = FALSE,
                             maf_data,
                             seq_type = "genome",
                             projection = "grch37",
                             from_indexed_flatfile = TRUE,
                             augmented = TRUE,
                             min_read_support = 3,
                             mode = "slms-3"){

  tabix_bin = "/home/rmorin/miniconda3/bin/tabix"
  table_name = config::get("results_tables")$ssm
  db = config::get("database_name")

  if(from_indexed_flatfile){
    base_path = config::get("project_base")

    #test if we have permissions for the full gambl + icgc merge
    if(mode == "slms-3"){
      if(augmented){
        maf_partial_path = config::get("results_flatfiles")$ssm$template$merged$augmented
      }else{
        maf_partial_path = config::get("results_flatfiles")$ssm$template$merged$deblacklisted
      }
    }else if (mode == "strelka2"){
      maf_partial_path = config::get("results_flatfiles")$ssm$all$strelka2
    }else{
      stop("You requested results from indexed flatfile. The mode should be set to either slms-3 (default) or strelka2. Please specify one of these modes.")
    }

    #setup file paths
    maf_path = glue::glue(maf_partial_path)
    full_maf_path = paste0(base_path, maf_path)
    full_maf_path_comp = paste0(base_path, maf_path, ".bgz")

  }

  if(!region == ""){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))

    if(projection == "grch37"){
      region = stringr::str_replace(region, "chr", "")
    }

    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = as.numeric(startend[1])
    qend = as.numeric(startend[2])
  }

  if(projection =="grch37"){
    chromosome = gsub("chr", "", chromosome)
  }

  if(missing(maf_data)){
    if(from_indexed_flatfile){
      #get column names for maf
      maf_head = as.vector(as.matrix(read.table(file = full_maf_path, header = FALSE, stringsAsFactors = FALSE, nrows = 1)))
      muts = system(paste(tabix_bin, full_maf_path_comp, region), intern = TRUE)

      if(length(muts) > 1){
        muts_region = vroom(I(muts), col_names = maf_head)
      }

      if(augmented){
        # drop poorly supported reads but only from augmented MAF
        muts_region = dplyr::filter(muts_region, t_alt_count >= min_read_support)
      }

    }else{
      con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
      muts_region = dplyr::tbl(con, table_name) %>%
        dplyr::filter(Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)

      muts_region = as.data.frame(muts_region)
      DBI::dbDisconnect(con)
    }
  }else{
    message("not using the database")
    muts_region = dplyr::filter(maf_data, Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)
    muts_region = dplyr::filter(maf_data, Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)
  }

  if(streamlined){
    muts_region = muts_region %>%
      dplyr::select(Start_Position, Tumor_Sample_Barcode)
  }else if(basic_columns){
    muts_region = muts_region[, c(1:45)]
  }

  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    muts_region = dplyr::select(muts_region, all_of(maf_cols))
  }

  #print all available columns
  if(!basic_columns && return_cols){
    print(colnames(muts_region))
  }

  muts_region$Chromosome = as.character(muts_region$Chromosome)

  return(muts_region)
}




#' Retrieve all coding SSMs from the GAMBL database in MAF-like format.
#'
#' @param limit_cohort Supply this to restrict mutations to one or more cohorts in a list.
#' @param exclude_cohort  Supply this to exclude mutations from one or more cohorts in a list.
#' @param limit_pathology Supply this to restrict mutations to one pathology.
#' @param limit_samples Supply this to restrict mutations to a vector of sample_id (instead of subsetting using the provided metadata)
#' @param these_samples_metadata Supply a metadata table to auto-subset the data to samples in that table before returning
#' @param force_unmatched_samples Optional argument for forcing unmatched samples, using get_ssm_by_samples.
#' @param projection Reference genome build for the coordinates in the MAF file. The default is hg19 genome build.
#' @param seq_type The seq_type you want back, default is genome.
#' @param basic_columns Set to TRUE to override the default behavior of returning only the first 45 columns of MAF data.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a list of indexes (integer) or a list of characters (matching columns in MAF).
#' @param return_cols If set to TRUE, a vector with all avaialble column names will be returned. Default is FALSE.
#' @param from_flatfile Set to TRUE to obtain mutations from a local flatfile instead of the database. This can be more efficient and is currently the only option for users who do not have ICGC data access.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs)
#' @param groups Unix groups for the samples to be included. Default is both gambl and icgc_dart samples.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is TRUE.
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' #basic usage
#' maf_data = get_coding_ssm(limit_cohort = c("BL_ICGC"))
#' maf_data = get_coding_ssm(limit_samples = "HTMCP-01-06-00485-01A-01D")
#'
get_coding_ssm = function(limit_cohort,
                          exclude_cohort,
                          limit_pathology,
                          limit_samples,
                          these_samples_metadata,
                          force_unmatched_samples,
                          projection = "grch37",
                          seq_type = "genome",
                          basic_columns = TRUE,
                          maf_cols = NULL,
                          return_cols = FALSE,
                          from_flatfile = TRUE,
                          augmented = TRUE,
                          min_read_support = 3,
                          groups = c("gambl", "icgc_dart"),
                          include_silent = TRUE){

  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  if(!missing(these_samples_metadata)){
    all_meta = these_samples_metadata
  }else{
    all_meta = get_gambl_metadata(from_flatfile = from_flatfile, seq_type_filter = seq_type)
  }

  all_meta = dplyr::filter(all_meta, seq_type == {{ seq_type }})

  #do all remaining filtering on the metadata then add the remaining sample_id to the query
  #unix groups
  all_meta = all_meta %>%
    dplyr::filter(unix_group %in% groups)

  #lmit cohort
  if(!missing(limit_cohort)){
    all_meta = all_meta %>%
      dplyr::filter(cohort %in% limit_cohort)
  }
  #exclude cohort
  if(!missing(exclude_cohort)){
    all_meta = all_meta %>%
      dplyr::filter(!cohort %in% exclude_cohort)
  }
  #limit pathology
  if(!missing(limit_pathology)){
    all_meta = all_meta %>%
      dplyr::filter(pathology %in% limit_pathology)
  }
  #limit samples
  if(!missing(limit_samples)){
    all_meta = all_meta %>%
      dplyr::filter(sample_id %in% limit_samples)
  }
  #pull info for loading .CDS.maf
  sample_ids = pull(all_meta, sample_id)

  #get file path for non-augmented maf
  if(from_flatfile && !augmented){
    maf_template = config::get("results_flatfiles")$ssm$template$cds$deblacklisted
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(config::get("project_base"), maf_path)
  }

  #get file path for augmented maf
  if(from_flatfile && augmented){
    maf_template = config::get("results_flatfiles")$ssm$template$cds$augmented
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(config::get("project_base"), maf_path)
  }

  #read file
  if(from_flatfile){
    message(paste("reading from:", full_maf_path))
    muts = fread_maf(full_maf_path) %>%
      dplyr::filter(Variant_Classification %in% coding_class) %>%
      as.data.frame()

    mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
    message(paste("mutations from", mutated_samples, "samples"))
  }

  #use db if not using flatfile
  if(!from_flatfile){
    table_name = config::get("results_tables")$ssm
    db = config::get("database_name")
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

    muts = tbl(con, table_name) %>%
      dplyr::filter(Variant_Classification %in% coding_class) %>%
      as.data.frame()

    DBI::dbDisconnect(con)
  }

  #if augmented maf selected, drop variants with low read support (default is 3)
  if(augmented){
    muts = dplyr::filter(muts, t_alt_count >= min_read_support)
  }

  #filter maf on selected sample ids
  muts = muts %>%
    dplyr::filter(Tumor_Sample_Barcode %in% sample_ids)

  mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
  message(paste("after linking with metadata, we have mutations from", mutated_samples, "samples"))

  #subset to fewer columns
  if(basic_columns){
    muts = muts[,c(1:45)]
  }

  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    muts = dplyr::select(muts, all_of(maf_cols))
  }

  #print all avaialble columns
  if(!basic_columns && return_cols){
    print(colnames(muts))
  }

  #drop rows for these samples so we can swap in the force_unmatched outputs instead
  if(!missing(force_unmatched_samples)){
    muts = muts %>%
      dplyr::filter(!sample_id %in% force_unmatched_samples)

    nsamp = length(force_unmatched_samples)
    message(paste("dropping variants from", nsamp, "samples and replacing with force_unmatched outputs"))

    #get replacements using get_ssm_by_samples
    fu_muts = get_ssm_by_samples(these_sample_ids = force_unmatched_samples)
    muts = bind_rows(muts, fu_muts)
  }

  return(muts)
}


#' Get the copy number and expression for a single gene.
#'
#' @param hugo_symbol One or more gene symbols. Should match the values in a maf file.
#' @param ensembl_gene_id One or more ensembl gene IDs. Only one of hugo_symbols or ensembl_gene_ids may be used.
#'
#' @return A list.
#' @import tidyverse
#' @export
#'
#' @examples
#' MYC_cn_expression = get_gene_cn_and_expression("MYC")
#'
get_gene_cn_and_expression = function(gene_symbol,
                                      ensembl_id){

    if(!missing(gene_symbol)){
      this_row = grch37_gene_coordinates %>%
        dplyr::filter(hugo_symbol == gene_symbol)

      this_region = paste0(this_row$chromosome, ":", this_row$start, "-", this_row$end)
      gene_name = gene_symbol
      }

    else{
      this_row = grch37_gene_coordinates %>%
        dplyr::filter(ensembl_gene_id == ensembl_id)

      this_region = paste0(this_row$chromosome, ":", this_row$start, "-",this_row$end)
      gene_name = ensembl_id
      gene_symbol = pull(this_row, hugo_symbol)
    }
  gene_cn = get_cn_states(regions_list = c(this_region), region_names = c(gene_name)) %>%
    as.data.frame()

  colnames(gene_cn)[1] = paste(colnames(gene_cn)[1], "CN", sep = "_")
  gene_cn = gene_cn %>%
    rownames_to_column("sample_id")

  gene_exp = get_gene_expression(hugo_symbols = c(gene_symbol), join_with = "genome")
  exp_copy = left_join(gene_cn, gene_exp, by = "sample_id")
  all_meta = get_gambl_metadata()
  exp_copy_meta = left_join(all_meta, exp_copy, by = "sample_id")
  return(exp_copy_meta)
}


#' Get the expression for one or more genes for all GAMBL samples
#'
#' @param metadata GAMBL metadata.
#' @param hugo_symbols One or more gene symbols. Should match the values in a maf file.
#' @param ensembl_gene_ids One or more ensembl gene IDs. Only one of hugo_symbols or ensembl_gene_ids may be used.
#' @param join_with How to restrict cases for the join. Can be one of genome, mrna or "any"
#' @param all_genes Set to TRUE to return full expression data frame without any subsetting. Avoid this if you don't want to use tons of RAM.
#' @param expression_data Optional argument to use an already loaded expression data frame (prevent function to re-load full df from flat file or database)
#' @param from_flatfile Deprecated but left here for backwards compatibility.
#'
#' @return A list.
#' @export
#'
#' @examples
#' MYC_expr = get_gene_expression(hugo_symbols = c("MYC"), join_with = "mrna")
#' Read full expression values df (no subsetting on genes)
#' full_expression_df = get_gene_expression_new(all_genes = TRUE, join_with = "genome")
#' Use loaded df (in previous step) to get expression values for IRF4 and MYC.
#' irf4_myc_expressions = get_gene_expression_new(hugo_symbols = c("IRF4", "MYC"), all_genes = FALSE, join_with = "genome", from_flatfile = FALSE, expression_data = full_expression_df)

get_gene_expression = function(metadata,
                               hugo_symbols,
                               ensembl_gene_ids,
                               join_with = "mrna",
                               all_genes = FALSE,
                               expression_data,
                               from_flatfile = TRUE){

  database_name = config::get("database_name")
  if(missing(metadata)){
    if(join_with == "mrna"){
      metadata = get_gambl_metadata(seq_type_filter = "mrna", only_available = FALSE)
      metadata = metadata %>%
        dplyr::select(sample_id)

      }else if(join_with == "genome"){
      metadata = get_gambl_metadata(only_available = FALSE)
      metadata = metadata %>%
        dplyr::select(sample_id)

      }else{
      metadata = get_gambl_metadata(seq_type_filter = c("genome","mrna"), only_available = FALSE)
      metadata = metadata %>%
        dplyr::select(sample_id, biopsy_id)
    }
  }

  if(missing(hugo_symbols) & missing(ensembl_gene_ids) & !all_genes){
    stop("ERROR: supply at least one gene symbol or Ensembl gene ID")
  }else if(!missing(hugo_symbols) & !missing(ensembl_gene_ids)){
    stop("ERROR: Both hugo_symbols and ensembl_gene_ids were provided. Please provide only one type of ID.")
  }
  tidy_expression_file = config::get("results_merged")$tidy_expression_file


  if(!missing(expression_data)){
    tidy_expression_data = as.data.frame(expression_data)
    if(!missing(hugo_symbols)){
      #lazily filter on the fly to conserve RAM
      wide_expression_data = tidy_expression_data %>%
        dplyr::filter(Hugo_Symbol %in% hugo_symbols) %>%
        dplyr::select(-ensembl_gene_id) %>%
        group_by(mrna_sample_id,Hugo_Symbol) %>% #deal with non 1:1 mapping of Hugo to Ensembl
        slice_head() %>%
        as.data.frame() %>%
        pivot_wider(names_from = Hugo_Symbol, values_from = expression)
    }
    if(!missing(ensembl_gene_ids)){
      wide_expression_data = tidy_expression_data %>%
        dplyr::filter(ensembl_gene_id %in% ensembl_gene_ids) %>%
        dplyr::select(-Hugo_Symbol) %>%
        as.data.frame() %>%
        pivot_wider(names_from = ensembl_gene_id, values_from = expression)
    }
  }else{
    #only ever load the full data frame when absolutely necessary
    if(all_genes & missing(ensembl_gene_ids) & missing(hugo_symbols)){
      tidy_expression_data = vroom::vroom(tidy_expression_file) %>%
        as.data.frame()
    }else{
      if(!missing(hugo_symbols)){
        #lazily filter on the fly to conserve RAM
        wide_expression_data = vroom::vroom(tidy_expression_file) %>%
          dplyr::filter(Hugo_Symbol %in% hugo_symbols) %>%
          dplyr::select(-ensembl_gene_id) %>%
          group_by(mrna_sample_id,Hugo_Symbol) %>% #deal with non 1:1 mapping of Hugo to Ensembl
          slice_head() %>%
          as.data.frame() %>%
          pivot_wider(names_from = Hugo_Symbol, values_from = expression)
      }
      if(!missing(ensembl_gene_ids)){
        wide_expression_data = vroom::vroom(tidy_expression_file) %>%
          dplyr::filter(ensembl_gene_id %in% ensembl_gene_ids) %>%
          dplyr::select(-Hugo_Symbol) %>%
          as.data.frame() %>%
          pivot_wider(names_from = ensembl_gene_id, values_from = expression)

      }

    }
  }

  if(join_with == "mrna" & missing(expression_data)){
    expression_wider = dplyr::select(wide_expression_data, -biopsy_id, -genome_sample_id)
    expression_wider = left_join(metadata, expression_wider, by = c("sample_id" = "mrna_sample_id"))

    }else if(join_with == "genome" & missing(expression_data)){
    expression_wider = dplyr::select(wide_expression_data, -mrna_sample_id, -biopsy_id) %>% dplyr::filter(genome_sample_id != "NA")
    expression_wider = left_join(metadata, expression_wider, by = c("sample_id" = "genome_sample_id"))

    }else if(join_with == "any" & missing(expression_data)){
    expression_wider = dplyr::select(wide_expression_data, -mrna_sample_id, -genome_sample_id)
    expression_wider = left_join(metadata, expression_wider, by = c("biopsy_id" = "biopsy_id"))

    }else if(join_with == "mrna" & !missing(expression_data)){
      expression_wider = wide_expression_data

    }else if(join_with == "genome" & !missing(expression_data)){
      expression_wider = wide_expression_data

    }else if(join_with == "any" & !missing(expression_data)){
      expression_wider = wide_expression_data
  }
  return(expression_wider)
}
