#global variables
coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")
cnames = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "NAME", "SOMATIC_SCORE", "STRAND_A", "STRAND_B", "TYPE", "FILTER", "VAF_tumour", "VAF_normal", "DP_tumour", "DP_normal", "tumour_sample_id", "normal_sample_id", "pair_status")
maf_header = c("Hugo_Symbol"=1,"Entrez_Gene_Id"=2,"Center"=3,"NCBI_Build"=4,"Chromosome"=5,"Start_Position"=6,"End_Position"=7,"Strand"=8,"Variant_Classification"=9,"Variant_Type"=10,"Reference_Allele"=11,"Tumor_Seq_Allele1"=12,"Tumor_Seq_Allele2"=13,"dbSNP_RS"=14,"dbSNP_Val_Status"=15,"Tumor_Sample_Barcode"=16,"Matched_Norm_Sample_Barcode"=17,"Match_Norm_Seq_Allele1"=18,"Match_Norm_Seq_Allele2"=19,"Tumor_Validation_Allele1"=20,"Tumor_Validation_Allele2"=21,"Match_Norm_Validation_Allele1"=22,"Match_Norm_Validation_Allele2"=23,"Verification_Status"=24,"Validation_Status"=25,"Mutation_Status"=26,"Sequencing_Phase"=27,"Sequence_Source"=28,"Validation_Method"=29,"Score"=30,"BAM_File"=31,"Sequencer"=32,"Tumor_Sample_UUID"=33,"Matched_Norm_Sample_UUID"=34,"HGVSc"=35,"HGVSp"=36,"HGVSp_Short"=37,"Transcript_ID"=38,"Exon_Number"=39,"t_depth"=40,"t_ref_count"=41,"t_alt_count"=42,"n_depth"=43,"n_ref_count"=44,"n_alt_count"=45,"all_effects"=46,"Allele"=47,"Gene"=48,"Feature"=49,"Feature_type"=50,"Consequence"=51,"cDNA_position"=52,"CDS_position"=53,"Protein_position"=54,"Amino_acids"=55,"Codons"=56,"Existing_variation"=57,"ALLELE_NUM"=58,"DISTANCE"=59,"STRAND_VEP"=60,"SYMBOL"=61,"SYMBOL_SOURCE"=62,"HGNC_ID"=63,"BIOTYPE"=64,"CANONICAL"=65,"CCDS"=66,"ENSP"=67,"SWISSPROT"=68,"TREMBL"=69,"UNIPARC"=70,"RefSeq"=71,"SIFT"=72,"PolyPhen"=73,"EXON"=74,"INTRON"=75,"DOMAINS"=76,"AF"=77,"AFR_AF"=78,"AMR_AF"=79,"ASN_AF"=80,"EAS_AF"=81,"EUR_AF"=82,"SAS_AF"=83,"AA_AF"=84,"EA_AF"=85,"CLIN_SIG"=86,"SOMATIC"=87,"PUBMED"=88,"MOTIF_NAME"=89,"MOTIF_POS"=90,"HIGH_INF_POS"=91,"MOTIF_SCORE_CHANGE"=92,"IMPACT"=93,"PICK"=94,"VARIANT_CLASS"=95,"TSL"=96,"HGVS_OFFSET"=97,"PHENO"=98,"MINIMISED"=99,"GENE_PHENO"=100,"FILTER"=101,"flanking_bps"=102,"vcf_id"=103,"vcf_qual"=104,"gnomAD_AF"=105,"gnomAD_AFR_AF"=106,"gnomAD_AMR_AF"=107,"gnomAD_ASJ_AF"=108,"gnomAD_EAS_AF"=109,"gnomAD_FIN_AF"=110,"gnomAD_NFE_AF"=111,"gnomAD_OTH_AF"=112,"gnomAD_SAS_AF"=113,"vcf_pos"=114,"gnomADg_AF"=115,"blacklist_count"=116)


#' @title Get Excluded Samples. 
#'
#' @description Exclude samples that have been excluded from certain analyses and drop from merges.
#'
#' @details Specify the tool or pipeline responsible for generating the files with `tool_name` and this function will return a vector of excluded sample IDs.
#'
#' @param tool_name The tool or pipeline that generated the files (should be the same for all). Default is "slms-3".
#'
#' @return A vector of sample IDs.
#'
#' @import readr dplyr
#' @export
#'
#' @examples
#' excluded_samp = get_excluded_samples()
#'
get_excluded_samples = function(tool_name = "slms-3"){
  base = check_config_value(config::get("repo_base"))

  #check for missingness
  path = paste0(base,"config/exclude.tsv")
  if(!file.exists(path)){
    print(paste("missing: ", path))
    message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
  }

  excluded_df = suppressMessages(read_tsv(paste0(base,"config/exclude.tsv")))
  excluded_samples = dplyr::filter(excluded_df, pipeline_exclude == tool_name) %>%
    pull(sample_id)

  return(excluded_samples)
}


#' @title Get SSM By Patients.
#'
#' @description Get MAF-format data frame for more than one patient.
#'
#' @details This function returns variants from a set of patients avoiding duplicated mutations from multiple samples from that patient (i.e. unique superset of variants).
#' This is done either by combining the contents of individual MAF files or subset from a merged MAF (wraps get_ssm_by_samples).
#' In most situations, this should never need to be run with `subset_from_merge = TRUE`. Instead use one of `get_coding_ssm` or `get_ssm_by_region`.
#' This function expects either a vector of patient IDs (`thse_patients_ids`) or an already subset metadata table (`these_samples_metadata`).
#'
#' @param these_patient_ids A vector of patient IDs that you want results for. The user can also use a metadata table that has been subset to the patient IDs of interest (`these_samples_metadata`).
#' @param these_samples_metadata A metadata subset to contain the rows corresponding to the patients of interest. If the vector of patient IDs is missing (`these_patient_ids`), this function will default to all patient IDs in the metadata table given to this parameter. 
#' @param tool_name Only supports slms-3 currently.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param seq_type The seq type you want results for. Default is "genome".
#' @param flavour Currently this function only supports one flavour option but this feature is meant for eventual compatibility with additional variant calling parameters/versions.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 43 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters (matching columns in MAF).
#' @param subset_from_merge Instead of merging individual MAFs, the data will be subset from a pre-merged MAF of samples with the specified `seq_type`.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead.
#' @param engine Specify one of readr or fread_maf (default) to change how the large files are loaded prior to subsetting. You may have better performance with one or the other but for me fread_maf is faster and uses a lot less RAM.
#'
#' @return A data frame with SSM calls for the selected patients in MAF format.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' #example 1, using a vector of patient IDs.
#' patients = c("00-14595", "00-15201", "01-12047")
#' patients_maf = get_ssm_by_patients(these_patient_ids = patients, seq_type = "genome", subset_from_merge = FALSE)
#'
#' #example 2, using a metadata table, subset to the patient IDs of interest.
#' patient_meta = get_gambl_metadata(seq_type_filter = "genome") %>% dplyr::filter(patient_id %in% patients)
#' patients_maf_2 = get_ssm_by_patients(these_samples_metadata = patient_meta,subset_from_merge = FALSE)
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
                               subset_from_merge = FALSE,
                               augmented = TRUE,
                               engine='fread_maf'){

  check_remote_configuration(auto_connect = TRUE)
  if(!subset_from_merge){
    message("WARNING: on-the-fly merges can be extremely slow and consume a lot of memory. Use at your own risk. ")
  }
  if(missing(these_patient_ids)){
    if(missing(these_samples_metadata)){
      stop("must supply either a vector of patient_ids or the metadata for those patients as these_samples_metadata")
    }
    these_patient_ids = pull(these_samples_metadata,patient_id) %>% unique()
  }
  augmented = TRUE
  #always requires augmented MAFs to ensure all variants from the patient are included
  to_exclude = get_excluded_samples(tool_name)
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter = seq_type) %>%
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
                            basic_columns = basic_columns,
                            maf_cols = maf_cols,
                            engine=engine))
}

#' @title Get SSM By Samples.
#'
#' @description Get MAF-format data frame for more than one sample and combine them together.
#'
#' @details This function internally runs `get_ssm_by_sample`. 
#' The user can either give the function a vector of sample IDs of interest with `these_sample_ids`,
#' or use a metadata table (`these_samples_metadata`), already subset to the sample IDs of interest.
#' In most situations, this should never need to be run with subset_from_merge = TRUE. 
#' Instead use one of `get_coding_ssm` or `get_ssm_by_region`.
#' See `get_ssm_by_sample` for more information.
#'
#' @param these_sample_ids A vector of sample_id that you want results for.
#' @param these_samples_metadata Optional metadata table. If provided, the function will return SSM calls for the sample IDs in the provided metadata table. 
#' @param tool_name Only supports slms-3 currently.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param seq_type  The seq type you want results for. Default is "genome".
#' @param flavour Currently this function only supports one flavour option but this feature is meant for eventual compatibility with additional variant calling parameters/versions.
#' @param these_genes A vector of genes to subset ssm to.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 45 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters.
#' @param subset_from_merge Instead of merging individual MAFs, the data will be subset from a pre-merged MAF of samples with the specified seq_type.
#' @param engine Specify one of readr or fread_maf (default) to change how the large files are loaded prior to subsetting. You may have better performance with one or the other but for me fread_maf is faster and uses a lot less RAM.
#'
#' @return A data frame in MAF format.
#'
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' #examples using the these_sample_ids parameter.
#' sample_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"))
#' hg38_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"),projection="hg38")
#' readr_sample_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"),subset_from_merge=TRUE,engine="readr")
#' slow_sample_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"),subset_from_merge=TRUE)
#'
#' #example using a metadata table subset to sample IDs of interest.
#' my_metadata = get_gambl_metadata(seq_type_filter = "genome") %>%
#'  dplyr::filter(pathology == "FL")
#'
#' sample_ssms = get_ssm_by_samples(these_samples_metadata = my_metadata)
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
                              subset_from_merge = FALSE,
                              augmented = TRUE,
                              engine = 'fread_maf'){

  remote_session = check_remote_configuration(auto_connect = TRUE)
  if(!subset_from_merge){
    message("WARNING: on-the-fly merges can be extremely slow and consume a lot of memory if many samples are involved. Use at your own risk. ")
  }
  to_exclude = get_excluded_samples(tool_name)

  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter = seq_type) %>%
      dplyr::filter(sample_id %in% these_sample_ids) %>%
      dplyr::filter(!sample_id %in% to_exclude)
  }else{
    if(missing(these_sample_ids)){
      #assume the user just wants the data for all the sample ids in this data frame
      these_sample_ids = pull(these_samples_metadata,sample_id)
    }
    else{
      these_samples_metadata = these_samples_metadata %>%
        dplyr::filter(sample_id %in% these_sample_ids) %>%
        dplyr::filter(!sample_id %in% to_exclude)
    }
  }
  #ensure we only have sample_id that are in the remaining metadata (no excluded/unavailable samples)
  these_sample_ids = these_sample_ids[which(these_sample_ids %in% these_samples_metadata$sample_id)]
  maf_column_types = "ccccciiccccccccccccccccccccccnccccccccciiiiii" #for the first 45 standard columns
  if(flavour=="legacy"){
    warning("I lied. Access to the old variant calls is not currently supported in this function")
    # TODO: implement loading of the old merged MAF under icgc_dart... vcf2maf-1.2 ..level_3 as per the other from_flatfile functions
    return()

  }else if(flavour=="clustered"){
    if(subset_from_merge && !augmented){
      maf_template = check_config_value(config::get("results_flatfiles")$ssm$template$merged$deblacklisted)
      maf_path = glue::glue(maf_template)
      full_maf_path =  paste0(check_config_value(config::get("project_base")), maf_path)
      message(paste("using existing merge:", full_maf_path))

      #check for missingness
      if(!file.exists(full_maf_path)){
        print(paste("missing: ", full_maf_path))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
        }

      if(engine=="fread_maf"){
        if(basic_columns){
          maf_df_merge = fread_maf(full_maf_path,select_cols = c(1:45)) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }else{
          maf_df_merge = fread_maf(full_maf_path) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }
      }else if(engine=="readr"){
        if(basic_columns){
          maf_df_merge = suppressMessages(read_tsv(full_maf_path,col_select = c(1:45),num_threads=12,col_types = maf_column_types,lazy = TRUE)) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }else{
          maf_df_merge = fread_maf(full_maf_path) %>%
            dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
            dplyr::filter(t_alt_count >= min_read_support)
        }
      }else{
        stop("specify one of readr or fread_maf as the file-reading engine")
      }

      if(!is.null(maf_cols) && !basic_columns){
        maf_df_merge = dplyr::select(maf_df_merge, all_of(maf_cols))
      }
    }

    if(subset_from_merge && augmented){
      maf_template = check_config_value(config::get("results_flatfiles")$ssm$template$merged$augmented)
      maf_path = glue::glue(maf_template)
      full_maf_path =  paste0(check_config_value(config::get("project_base")), maf_path)
      message(paste("using existing merge:", full_maf_path))

      #check for missingness
      if(!file.exists(full_maf_path)){
        print(paste("missing: ", full_maf_path))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
      }

      #maf_df_merge = read_tsv(full_maf_path) %>%
      #  dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
      #  dplyr::filter(t_alt_count >= min_read_support)
      if(basic_columns){
        maf_df_merge = fread_maf(full_maf_path,select_cols = c(1:45)) %>%
          dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
          dplyr::filter(t_alt_count >= min_read_support)
      }else{
        maf_df_merge = fread_maf(full_maf_path) %>%
          dplyr::filter(Tumor_Sample_Barcode %in% these_sample_ids) %>%
          dplyr::filter(t_alt_count >= min_read_support)
      }
      #subset maf to only include first 43 columns (default)
      if(basic_columns){maf_df_merge = dplyr::select(maf_df_merge, c(1:45))}
      #subset maf to a specific set of columns (defined in maf_cols)
      if(!is.null(maf_cols) && !basic_columns){maf_df_merge = dplyr::select(maf_df_merge, all_of(maf_cols))}
    }

    if(!subset_from_merge){
      if(remote_session){
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
            verbose = FALSE)
          maf_df_list[[this_sample]]=maf_df
        }
      }else{
        maf_df_list = mclapply(these_sample_ids,function(x){get_ssm_by_sample(
        this_sample_id=x,
        these_samples_metadata = these_samples_metadata,
        tool_name = tool_name,
        projection = projection,
        augmented = augmented,
        flavour = flavour,
        min_read_support = min_read_support,
        basic_columns = basic_columns,
        maf_cols = maf_cols,
        verbose = FALSE
        )},mc.cores = 12)
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


#' @title Get SSM By Sample.
#'
#' @description Get the SSMs (i.e. load MAF) for a single sample.
#'
#' @details This was implemented to allow flexibility because there are some samples that we may want to use a different set of variants than those in the main GAMBL merge.
#' The current use case is to allow a force_unmatched output to be used to replace the SSMs from the merge for samples with known contamination in the normal.
#' This will also be useful to apply a blacklist to individual MAFs when coupled with annotate_ssm_blacklist.
#'
#' @param this_sample_id Required. The sample_id you want the data from.
#' @param this_seq_type Required if not specifying these_samples_metadata. The seq_type of the sample you want data from.
#' @param these_samples_metadata Required if not specifying both this_sample_id and this_seq_type a single row or entire metadata table containing your sample_id.
#' @param tool_name The name of the variant calling pipeline (currently only slms-3 is supported).
#' @param projection The projection genome build. Supports hg38 and grch37.
#' @param these_genes A vector of genes to subset ssm to.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF.
#' @param flavour Currently this function only supports one flavour option but this feature is meant for eventual compatibility with additional variant calling parameters/versions.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Return first 43 columns of MAF rather than full details. Default is TRUE.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters.
#' @param verbose Enable for debugging/noisier output.
#'
#' @return data frame in MAF format.
#'
#' @import dplyr tidyr
#' @export
#'
#' @examples
#' this_sample_df = get_ssm_by_sample(this_sample_id = "HTMCP-01-06-00485-01A-01D", this_seq_type = "genome",tool_name = "slims-3", projection = "grch37")
#' capture_meta = get_gambl_metadata(seq_type_filter = "capture")
#' ssm_sample = get_ssm_by_sample(this_sample_id = "CASA0002_2015-03-10", projection = "grch37",augmented = T,these_samples_metadata = capture_meta)
#'
get_ssm_by_sample = function(this_sample_id,
                             this_seq_type,
                             these_samples_metadata,
                             tool_name = "slms-3",
                             projection = "grch37",
                             these_genes,
                             augmented = TRUE,
                             flavour = "clustered",
                             min_read_support = 3,
                             basic_columns = TRUE,
                             maf_cols = NULL,
                             verbose = FALSE
                             ){
  remote_session = check_remote_configuration(auto_connect = TRUE)
  if(missing(this_seq_type) & missing(these_samples_metadata)){
    stop("Must provide both a sample_id and seq_type for that sample via this_sample_id and this_seq_type")
  }
  if(missing(this_seq_type)){
    #get it from the metadata
    this_seq_type = dplyr::filter(these_samples_metadata,sample_id==this_sample_id) %>% pull(seq_type)
  }
  #figure out which unix_group this sample belongs to
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter = this_seq_type) %>%
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
    path_template = check_config_value(config::get("results_flatfiles",config="default")$ssm$template$clustered$deblacklisted)
    path_complete = unname(unlist(glue::glue(path_template)))
    full_maf_path = paste0(check_config_value(config::get("project_base",config="default")), path_complete)
    local_full_maf_path = paste0(check_config_value(config::get("project_base")), path_complete)
    if(augmented){
      path_template = check_config_value(config::get("results_flatfiles",config="default")$ssm$template$clustered$augmented)
      path_complete = unname(unlist(glue::glue(path_template)))
      aug_maf_path = paste0(check_config_value(config::get("project_base",config="default")), path_complete)
      local_aug_maf_path = paste0(check_config_value(config::get("project_base")), path_complete)
    }
  }else{
    warning("Currently the only flavour available to this function is 'clustered'")
  }
  if(remote_session){
    #check if file exists
    status = ssh::ssh_exec_internal(ssh_session,command=paste("stat",aug_maf_path),error=F)$status
    #aug_maf_path = paste0(aug_maf_path,".gz")
    #local_aug_maf_path = paste0(local_aug_maf_path,".gz")
    #full_maf_path = paste0(full_maf_path,".gz")
    #local_full_maf_path = paste0(local_full_maf_path,".gz")
    #deprecate the usage of gzipped MAF for now

    # first check if we already have a local copy
    # Load data from local copy or get a local copy from the remote path first
    if(status==0){
      if(verbose){
        print(paste("found:",aug_maf_path))
        print(paste("local home:",local_aug_maf_path))
      }
      dirN = dirname(local_aug_maf_path)

      suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
      if(!file.exists(local_aug_maf_path)){

        ssh::scp_download(ssh_session,aug_maf_path,dirN)
      }

     #check for missingness
     if(!file.exists(local_aug_maf_path)){
      print(paste("missing: ", local_aug_maf_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
     }

      sample_ssm = fread_maf(local_aug_maf_path) %>%
      dplyr::filter(t_alt_count >= min_read_support)
    }else{
      if(verbose){
        print(paste("will use",full_maf_path))
        print(paste("local home:",local_full_maf_path))
      }
      dirN = dirname(local_full_maf_path)

      suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
      if(!file.exists(local_full_maf_path)){

        ssh::scp_download(ssh_session,full_maf_path,dirN)
      }
      #check for missingness
      if(!file.exists(local_full_maf_path)){
        print(paste("missing: ", local_full_maf_path))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote"')
      }

      sample_ssm = fread_maf(local_full_maf_path)
    }
  }else if(augmented && file.exists(aug_maf_path)){
    full_maf_path = aug_maf_path

    #check for missingness
    if(!file.exists(full_maf_path)){
      print(paste("missing: ",full_maf_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    sample_ssm = fread_maf(full_maf_path)
    if(min_read_support){
      # drop poorly supported reads but only from augmented MAF
      sample_ssm = dplyr::filter(sample_ssm, t_alt_count >= min_read_support)
    }
  }else{
    if(!file.exists(full_maf_path)){
      print(paste("missing: ", full_maf_path))
      message(paste("warning: file does not exist, skipping it.", full_maf_path))
      return()
    }
    #check for missingness
    if(!file.exists(full_maf_path)){
      print(paste("missing: ", full_maf_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
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

  return(sample_ssm)
}


#' @title Get GAMBL metadata.
#'
#' @description Return metadata for a selection of samples.
#'
#' @details This function returns metadata for GAMBL samples. Options for subset and filter the returned data are available.
#' For more information on how to use this function with different filtering criteria, refer to the parameter descriptions,
#' examples and vignettes. Embargoed cases (current options: 'BLGSP-study', 'FL-study', 'DLBCL-study', 'FL-DLBCL-study', 'FL-DLBCL-all', 'DLBCL-unembargoed', 'BL-DLBCL-manuscript', 'MCL','MCL-CLL')
#' 
#' @param seq_type_filter Filtering criteria (default: all genomes).
#' @param tissue_status_filter Filtering criteria (default: only tumour genomes, can be "mrna" or "any" for the superset of cases).
#' @param case_set Optional short name for a pre-defined set of cases avoiding any embargoed cases (current options: 'BLGSP-study', 'FL-study', 'DLBCL-study', 'FL-DLBCL-study', 'FL-DLBCL-all', 'DLBCL-unembargoed', 'BL-DLBCL-manuscript', 'MCL','MCL-CLL').
#' @param remove_benchmarking By default the FFPE benchmarking duplicate samples will be dropped.
#' @param sample_flatfile Optionally provide the full path to a sample table to use instead of the default.
#' @param biopsy_flatfile Optionally provide the full path to a biopsy table to use instead of the default.
#' @param with_outcomes Optionally join to gambl outcome data.
#' @param only_available If TRUE, will remove samples with FALSE or NA in the bam_available column (default: TRUE).
#' @param from_flatfile New default is to use the metadata in the flat-files from your clone of the repo. Can be overridden to use the database.
#' @param seq_type_priority For duplicate sample_id with different seq_type available, the metadata will prioritize this seq_type and drop the others.
#'
#' @return A data frame with metadata for each biopsy in GAMBL
#'
#' @import config dplyr tidyr readr RMariaDB DBI
#' @export
#'
#' @examples
#' #basic usage 
#' my_metadata = get_gambl_metadata()
#'
#' #use pre-defined custom sample sets
#' only_blgsp_metadata = get_gambl_metadata(case_set="BLGSP-study")
#'
#' #override default filters and request metadata for samples other than tumour genomes, e.g. also get the normals
#' only_normal_metadata = get_gambl_metadata(tissue_status_filter = c('tumour','normal'))
#'
#' non_duplicated_genome_and_capture = get_gambl_metadata(seq_type_filter=c('genome','capture'),seq_type_priority="genome")
#'
get_gambl_metadata = function(seq_type_filter = "genome",
                              tissue_status_filter = "tumour",
                              case_set,
                              remove_benchmarking = TRUE,
                              with_outcomes = TRUE,
                              from_flatfile = TRUE,
                              sample_flatfile,
                              biopsy_flatfile,
                              only_available = TRUE,
                              seq_type_priority = "genome"){

  check_remote_configuration()
  #this needs to be in any function that reads files from the bundled GAMBL outputs synced by Snakemake
  outcome_table = get_gambl_outcomes(from_flatfile = from_flatfile)

  if(from_flatfile){
    base = config::get("repo_base")
    if(missing(sample_flatfile)){
      sample_flatfile = paste0(base, config::get("table_flatfiles")$samples)
    }
    if(missing(biopsy_flatfile)){
      biopsy_flatfile = paste0(base, config::get("table_flatfiles")$biopsies)
    }

    #check for missingness
    if(!file.exists(sample_flatfile)){
      print(paste("missing: ", sample_flatfile))
      message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
    }

    if(!file.exists(biopsy_flatfile)){
      print(paste("missing: ", biopsy_flatfile))
      message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
    }

    sample_meta = suppressMessages(read_tsv(sample_flatfile, guess_max = 100000))
    biopsy_meta = suppressMessages(read_tsv(biopsy_flatfile, guess_max = 100000))

  }else{
    db = check_config_value(config::get("database_name"))
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
    dplyr::filter(seq_type %in% seq_type_filter & tissue_status %in% tissue_status_filter) %>%
    dplyr::select(-sex)
  
  #if only normals were requested, just return what we have because there is nothing else to join
  #if(tissue_status_filter == "normal"){
  #  return(sample_meta)
  #}

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
    case_set_path = check_config_value(config::get("sample_sets")$default)
    full_case_set_path =  paste0(check_config_value(config::get("repo_base")), case_set_path)
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
      transformed_manual <- paste0(
          base,
          "data/metadata/raw_metadata/gambl_tFL_manual.tsv"
      )
      transformed_manual <- suppressMessages(
          read_tsv(
              transformed_manual
          )
      )

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

      gambl_transformations <- paste0(
          base,
          "data/metadata/raw_metadata/gambl_transformation.txt"
      )
      gambl_transformations <- suppressMessages(
              read_delim(
                  gambl_transformations,
                  delim = " "
              )
          ) %>%
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

      curated <- paste0(
          base,
          "data/metadata/raw_metadata/clin_review_fl.tsv"
      )
      curated <- suppressMessages(
          read_tsv(
              curated
          )
      )

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
      set_file = paste0(base, "data/metadata/BLGSP--DLBCL-case-set.tsv")
      adult_bl_manuscript_samples = read_tsv(set_file) %>%
        pull(Tumor_Sample_Barcode)

      all_meta = all_meta %>%
      dplyr::filter(sample_id %in% adult_bl_manuscript_samples)

    }else if(case_set == "BL-DLBCL-manuscript-HTMCP"){
      set_file = paste0(base, "data/metadata/BLGSP--DLBCL-case-set.tsv")
      adult_bl_manuscript_samples = read_tsv(set_file) %>%
        pull(Tumor_Sample_Barcode)

      all_meta = all_meta %>%
      dplyr::filter(sample_id %in% adult_bl_manuscript_samples | cohort == "DLBCL_HTMCP")

    }else if(case_set == "FL-DLBCL-all"){
      set_file = paste0(base, "data/metadata/FL--DLBCL--all-case-set.tsv")
      fl_dlbcl_all_samples = read_tsv(set_file) %>%
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

  all_meta = GAMBLR::tidy_lymphgen(all_meta,
              lymphgen_column_in = "lymphgen_cnv_noA53",
              lymphgen_column_out = "lymphgen",
              relevel=TRUE)

  #all_meta = GAMBLR::collate_lymphgen(all_meta, verbose=FALSE)

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

#' @title Add ICGC metadata.
#'
#' @description Layer on ICGC metadata from a supplemental table to fill in missing COO.
#'
#' @details INTERNAL FUNCTION called by `get_gambl_metadata`, not meant for out-of-package usage.
#'
#' @param incoming_metadata A metadata table (probably output from `get_gambl_metadata`).
#'
#' @return Metadata with layered information (ICGC).
#'
#' @import dplyr readr stringr
#'
#' @examples
#' icgc_metadata = add_icgc_metadata(incoming_metadata = my_meta)
#'
add_icgc_metadata = function(incoming_metadata){
  repo_base = check_config_value(config::get("repo_base"))
  icgc_publ_file = paste0(repo_base,"data/metadata/raw_metadata/MALY_DE_tableS1.csv")
  icgc_publ = suppressMessages(suppressWarnings(read_csv(icgc_publ_file)))
  icgc_publ = icgc_publ[,c(1:20)]
  #fix commas as decimals
  icgc_publ = mutate(icgc_publ, purity = str_replace(purity, ",", "."))
  icgc_publ = mutate(icgc_publ, sex = str_to_upper(sex))

  icgc_raw_path = paste0(repo_base,"data/metadata/raw_metadata/ICGC_MALY_seq_md.tsv")

  #check for missingness
  if(!file.exists(icgc_raw_path)){
    print(paste("missing: ", icgc_raw_path))
    message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
  }

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


#' @title Get GAMBL Outcomes.
#'
#' @description Get the patient-centric clinical metadata.
#'
#' @details INTERNAL FUNCTION called by `get_gambl_metadata`, not meant for out-of-package usage.
#'
#' @param patient_ids Vector of patient IDs.
#' @param time_unit Return follow-up times in one of three time units: year, month or day. Default is "year".
#' @param censor_cbioportal Optionally request the censoring to be encoded in the specific style required by cBioPortal. Default is FALSE.
#' @param complete_missing Optionally fill in any gaps to ensure we have values for every patient (censor at 0 if missing). Default is FALSE.
#' @param from_flatfile Optionally set to FALSE to use the database to get the survival data. Default is TRUE.
#'
#' @return Data frame with one row for each patient_id.
#'
#' @import tidyr dplyr readr RMariaDB DBI
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
    outcome_flatfile = paste0(check_config_value(config::get("repo_base")), check_config_value(config::get("table_flatfiles")$outcomes))

    #check for missingness
    if(!file.exists(outcome_flatfile)){
      print(paste("missing: ", outcome_flatfile))
      message("Have you cloned the GAMBL repo and added the path to this directory under the local section of your config?")
    }

    all_outcome = suppressMessages(read_tsv(outcome_flatfile))

  }else{
    db = check_config_value(config::get("database_name"))
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


#' @title Get Combined SV.
#'
#' @description Retrieve Combined Manta and GRIDSS-derived SVs from a flatfile and filter.
#'
#' @details The bedpe files used as input to this function were pre-filtered for a minimum VAF of 0.05, and SVs affecting.
#' common translocation regions (BCL2, BCL6, MYC, CCND1) were whitelisted (e.g. no VAF filter applied).
#' Therefore if you wish to post-filter the SVs we recommend doing so carefully after loading this data frame.
#' Further, the input bedpe file is annotated with oncogenes and superenhancers from naive and germinal centre B-cells.
#' You can subset to events affecting certain loci using the "oncogenes" argument.
#'
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Recommended: 0. (default: 0)
#' @param these_sample_ids A character vector of tumour sample IDs you wish to retrieve SVs for.
#' @param with_chr_prefix Prepend all chromosome names with chr (required by some downstream analyses). Default is FALSE.
#' @param projection The projection genome build. Default is "grch37".
#' @param oncogenes A character vector of genes commonly involved in translocations. Possible values: CCND1, CIITA, SOCS1, BCL2, RFTN1, BCL6, MYC, PAX5.
#'
#' @return A data frame in a bedpe-like format with additional columns that allow filtering of high-confidence SVs.
#'
#' @import config dplyr readr stringr
#' @export
#'
#' @examples
#' get_combined_sv(oncogenes = c("MYC", "BCL2", "BCL6"))
#'
get_combined_sv = function(min_vaf = 0,
                           these_sample_ids,
                           with_chr_prefix = FALSE,
                           projection = "grch37",
                           oncogenes){

  base_path = check_config_value(config::get("project_base"))
  sv_file = check_config_value(config::get()$results_flatfiles$sv_combined$icgc_dart)
  if(projection == "hg38"){
    sv_file = str_replace(sv_file, "--grch37", "--hg38")
  }
  sv_file = paste0(base_path, sv_file)
  permissions = file.access(sv_file, 4)
  if(permissions == - 1){
    sv_file = check_config_value(config::get()$results_flatfiles$sv_combined$gambl)
    sv_file = paste0(base_path, sv_file)
  }

  #check for missingness
  if(!file.exists(sv_file)){
    print(paste("missing: ", sv_file))
    message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
    message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
  }

  all_sv = suppressMessages(read_tsv(sv_file, col_types = "cnncnncnccccnnccncn")) %>%
    dplyr::rename(c("VAF_tumour" = "VAF")) %>%
    dplyr::filter(VAF_tumour >= min_vaf)

  if(!missing(these_sample_ids)){
    all_sv = all_sv %>%
      dplyr::filter(tumour_sample_id %in% these_sample_ids)
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

  all_sv = all_sv %>%
      mutate(FILTER = "PASS") #i.e all variants returned with get_combined_sv() all have PASS in the FILTER column.

  return(all_sv)
}


#' @title Get Manta SVs
#'
#' @description Retrieve Manta SVs and filter.
#'
#' @details Return Manta SVs with aditional VCF information to allow for filtering of high-confidence variants.
#' To return SV calls for multiple samples, give `these_sample_ids` a vector of sample IDs, if only one sample is desired,
#' give this parameter one sample ID, as a string (or a vector of characters). The user can also call the `these_samples_metadata`
#' parameter to make use of an already subset metadata table. In this case, the returned calls will be restricted to the sample_ids 
#' within that data frame. This function relies on a set of specific functions to be successful in returning SV calls for any 
#' available sample in gambl. First, this function calls `get_combined_sv` and performs an `anit_join` with the full metadata to 
#' identify what samples are currently missing from the return of `get_combined_sv`. This function then calls `get_manta_sv_by_samples` 
#' (wrapper function for `get_manta_sv_by_sample`) on the subset of the missing samples. The merged calls are subject to any 
#' filtering that is specified within this function. This function can also restrict the returned calls to any genomic regions 
#' specified within `chromosome`, `qstart`, `qend`, or the complete region specified under `region` (in chr:start-end format). 
#' Useful filtering parameters are also available, use `min_vaf` to set the minimum tumour VAF for a SV to be returned and `min_score` 
#' to set the lowest Manta somatic score for a SV to be returned. `pair_status` can be used to only return variants that are 
#' annotated with PASS in the filtering column (VCF).
#'
#' @param these_sample_ids A vector of multiple sample_id (or a single sample ID as a string) that you want results for.
#' @param these_samples_metadata A metadata table to auto-subset the data to samples in that table before returning.
#' @param projection The projection genome build.
#' @param chromosome Optional, the chromosome you are restricting to.
#' @param qstart Optional, query start coordinate of the range you are restricting to.
#' @param qend Optional, query end coordinate of the range you are restricting to.
#' @param region Optional, region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Default is 0.1.
#' @param min_score The lowest Manta somatic score for a SV to be returned. Default is 40.
#' @param pass If set to TRUE, only return SVs that are annotated with PASS in the FILTER column. Set to FALSE to keep all variants, regardless if they PASS the filters. Default is TRUE. 
#' @param pairing_status Use to restrict results (if desired) to matched or unmatched results (default is to return all).
#' @param from_flatfile Set to TRUE by default, FALSE is no longer supported (database).
#'
#' @return A data frame in a bedpe-like format with additional columns that allow filtering of high-confidence SVs.
#' 
#' @import dplyr
#' @export
#'
#' @examples
#' #lazily get every SV in the table with default quality filters
#' all_sv = get_manta_sv()
#'
#' #get all SVs for a single sample
#' some_sv = get_manta_sv(these_sample_ids = "94-15772_tumorA")
#'
#' #get the SVs in a region around MYC
#' myc_locus_sv = get_manta_sv(region = "8:128723128-128774067")
#'
#' #get SVs for multiple samples, using these_samples_id
#' my_samples = get_gambl_metadata() %>% dplyr::select(sample_id) %>% head(10) %>% pull(sample_id)
#' my_svs_2 = get_manta_sv(these_sample_ids = my_samples, projection = "hg38")
#'
#' #get SVs for multiple samples using a metadata table and with no VAF/score filtering
#' my_metadata = get_gambl_metadata() %>% head(10)
#' my_svs = get_manta_sv(these_samples_metadata = my_metadata, min_vaf = 0, min_score = 0)
#'
get_manta_sv = function(these_sample_ids,
                        these_samples_metadata,
                        projection = "grch37",
                        chromosome,
                        qstart,
                        qend,
                        region,
                        min_vaf = 0.1,
                        min_score = 40,
                        pass = TRUE,
                        pairing_status,
                        from_flatfile = TRUE){
  
  if(!missing(region)){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }

  if(from_flatfile){
    all_sv = get_combined_sv(projection = projection)

    all_meta = get_gambl_metadata()

    #add pairing status to get_combined_sv return
    sub_meta = all_meta %>%
      dplyr::select(sample_id, pairing_status) %>%
      rename(pair_status = pairing_status)
      
    all_sv = left_join(all_sv, sub_meta, by = c("tumour_sample_id" = "sample_id"))

    #get metadata for samples currently missing from the merged results
    missing_samples = all_meta %>%
      anti_join(all_sv, by = c("sample_id" = "tumour_sample_id"))
    
    #call get manta_sv_by_samples on samples missing from current merge
    missing_sv = get_manta_sv_by_samples(these_samples_metadata = missing_samples,
                                         projection = projection,
                                         min_vaf = min_vaf,
                                         min_score = min_score,
                                         pass = pass)
    
    #combine current manta merged results with missing samples
    all_sv = bind_rows(all_sv, missing_sv)
    
  }else{
    stop("database usage is deprecated, please set from_flatfile to TRUE...")
  }
  
  if(!missing(region) || !missing(chromosome)){
    suppressWarnings({
      if(grepl("chr",chromosome)){
        chromosome = gsub("chr", "", chromosome)
      }
    })
    
    all_sv = all_sv %>%
      dplyr::filter((CHROM_A == chromosome & START_A >= qstart & START_A <= qend) | (CHROM_B == chromosome & START_B >= qstart & START_B <= qend))
  }
  
  #VAF and somatic score filtering
  all_sv = all_sv %>%
    dplyr::filter(VAF_tumour >= min_vaf & SCORE >= min_score)
  
  #PASS filter
  if(pass){
    all_sv = all_sv %>%
      dplyr::filter(FILTER == "PASS")
  }
  
  #pairing status filter
  if(!missing(pairing_status)){
    all_sv = all_sv %>%
      dplyr::filter(pair_status == pairing_status)
  }

  #sample IDs filter
  if(!missing(these_sample_ids) && missing(these_samples_metadata)){
    all_sv = all_sv %>%
      dplyr::filter(tumour_sample_id %in% these_sample_ids)
  }
  
  #metadata filter
  if(!missing(these_samples_metadata) && missing(these_sample_ids)){
    all_sv = all_sv %>%
      dplyr::filter(tumour_sample_id %in% these_samples_metadata$sample_id)
  }
  
  #as data frame
  all_sv = as.data.frame(all_sv)

  return(all_sv)
}


#' @title Get Lymphgen.
#'
#' @description Get a specific flavour of LymphGen from the main GAMBL outputs.
#'
#' @details Get a specific flavour of LymphGen from the main GAMBL outputs and tidy the composites.
#' Optionally return a matrix of features instead
#'
#' @param flavour Lymphgen flavour.
#' @param these_samples_metadata A metadata table to auto-subset the data to samples in that table before returning.
#' @param return_feature_matrix Boolean parameter, default is FALSE.
#' @param return_feature_annotation Boolean parameter, default is FALSE.
#' @param lymphgen_file Path to lymphgen file.
#' @param keep_all_rows Boolean parameter, default is FALSE.
#' @param keep_original_columns Boolean parameter, default is FALSE.
#'
#' @return A data frame.
#' 
#' @import config dplyr tidyr readr stringr tibble 
#' @export
#'
#' @examples
#' lymphgens = get_lymphgen(flavour = "no_cnvs.no_sv.with_A53")
#'
get_lymphgen = function(these_samples_metadata,
                        flavour,
                        return_feature_matrix = FALSE,
                        return_feature_annotation = FALSE,
                        lymphgen_file,
                        keep_all_rows = FALSE,
                        keep_original_columns = FALSE){

  if(missing(these_samples_metadata)){
    if(!keep_all_rows){
      these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome")
    }
  }
  if(missing(flavour)){
    if(!missing(lymphgen_file)){
      lg_path = lymphgen_file
    }else{
      message("please provide a path to your lymphgen output file or one of the following flavours")
      print(check_config_value(config::get("results_merged_wildcards")$lymphgen_template))
      return(NULL)
    }
  }else{
    lg_path = paste0(check_config_value(config::get("project_base")), check_config_value(config::get("results_merged")$lymphgen_template))
    lg_path = glue::glue(lg_path)
  }

  lg = suppressMessages(read_tsv(lg_path))
  lg_tidy = tidy_lymphgen(lg,lymphgen_column_in = "Subtype.Prediction",lymphgen_column_out = "LymphGen")
  if(return_feature_matrix | return_feature_annotation){
    lg_ord = select(lg_tidy,Sample.Name,LymphGen) %>% arrange(LymphGen) %>% pull(Sample.Name)
    lg_levels = select(lg_tidy,Sample.Name,LymphGen) %>% arrange(LymphGen) %>% pull(LymphGen)
    all_mcd = separate(lg_tidy,col="MCD.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>%
      pull(MCD) %>% unique()
    all_mcd_genes = str_remove(all_mcd,"_.*")%>% unique()
    all_mcd_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_mcd_genes)
    feat_mcd = separate(lg_tidy,col="MCD.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1)

    feat_mcd_genes = separate(lg_tidy,col="MCD.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()


    mcd_mat = left_join(all_mcd_df,feat_mcd_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_mcd = mutate(feat_mcd_genes,Class="MCD")

    all_ezb = separate(lg_tidy,col="EZB.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>% pull(MCD) %>% unique()
    all_ezb_genes = str_remove(all_ezb,"_.*")%>% unique()

    feat_ezb_genes = separate(lg_tidy,col="EZB.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()

    all_ezb_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_ezb_genes)
    feat_ezb = separate(lg_tidy,col="EZB.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1)

    ezb_mat = left_join(all_ezb_df,feat_ezb_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_ezb = mutate(feat_ezb_genes,Class="EZB")

    all_bn2 = separate(lg_tidy,col="BN2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>% pull(MCD) %>% unique()
    all_bn2_genes = str_remove(all_bn2,"_.*")%>% unique()
    all_bn2_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_bn2_genes)

    feat_bn2 = separate(lg_tidy,col="BN2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name,Feature) %>% mutate(present=1)

    feat_bn2_genes = separate(lg_tidy,col="BN2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()


    bn2_mat = left_join(all_bn2_df,feat_bn2_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_bn2 = mutate(feat_bn2_genes,Class="BN2")

    all_st2 = separate(lg_tidy,col="ST2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>%
      pull(MCD) %>% unique()
    all_st2_genes = str_remove(all_st2,"_.*") %>% unique()
    all_st2_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_st2_genes)

    feat_st2 = separate(lg_tidy,col="ST2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name,Feature) %>% mutate(present=1)

    feat_st2_genes = separate(lg_tidy,col="ST2.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()


    st2_mat = left_join(all_st2_df,feat_st2_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_st2 = mutate(feat_st2_genes,Class="ST2")

    all_n1 = separate(lg_tidy,col="N1.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "MCD") %>% dplyr::filter(!is.na(MCD)) %>% pull(MCD) %>% unique()
    all_n1_genes = str_remove(all_n1,"_.*") %>% unique()

    #all_n1_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_n1)
    all_n1_df = expand.grid(Sample.Name=unique(lg_tidy$Sample.Name),Feature=all_n1_genes)
    feat_n1 = separate(lg_tidy,col="N1.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>% dplyr::filter(!is.na(Feature)) %>%
      select(Sample.Name,Feature) %>% mutate(present=1)

    feat_n1_genes = separate(lg_tidy,col="N1.Features",into=c(paste0("Feature_MCD_",seq(1:15))),sep=",") %>%
      pivot_longer(starts_with("Feature_"),values_to = "Feature") %>%
      dplyr::filter(!is.na(Feature)) %>% select(Sample.Name,Feature) %>% mutate(present=1) %>%
      mutate(Feature=str_remove(Feature,"_.*")) %>% group_by(Sample.Name,Feature) %>% slice_head()

    #n1_mat = left_join(all_n1_df,feat_n1) %>% mutate(present=replace_na(present,0)) %>%
    #  pivot_wider(names_from="Feature",values_from="present")

    n1_mat = left_join(all_n1_df,feat_n1_genes) %>% mutate(present=replace_na(present,0)) %>%
      pivot_wider(names_from="Feature",values_from="present")
    feat_n1 = mutate(feat_n1_genes,Class="N1")

    all_genes = c(all_n1_genes,all_ezb_genes,all_st2_genes,all_bn2_genes,all_mcd_genes)
    print(table(all_genes))
    feat_all = bind_rows(feat_n1,feat_st2,feat_mcd,feat_ezb,feat_bn2)

  if(return_feature_annotation){
    #just give the user the association between each feature and its class along with some summary stats
    feat_count = group_by(feat_all,Feature,Class) %>% count()
    return(feat_count)
  }

  all_mat = left_join(ezb_mat,mcd_mat)
  all_mat = left_join(all_mat,bn2_mat)
  all_mat = left_join(all_mat,n1_mat)
  all_mat = left_join(all_mat,st2_mat)
  if(!keep_all_rows){
    all_mat = dplyr::filter(all_mat,Sample.Name %in% these_samples_metadata$sample_id)
  }
  all_mat = all_mat %>% column_to_rownames("Sample.Name")

    return(all_mat)
  }else{
    if(!keep_original_columns){
      lg_tidy = dplyr::select(lg_tidy,Sample.Name,LymphGen)
    }
    lg_tidy = lg_tidy %>% dplyr::rename("sample_id"="Sample.Name")
    if(!keep_all_rows){
      lg_tidy = dplyr::filter(lg_tidy,sample_id %in% these_samples_metadata$sample_id)
    }
    return(lg_tidy)

  }
}


#' @title Get CN States.
#'
#' @description Get a copy number matrix for all samples based on segmented data in the database.
#'
#' @details This function returns CN states for the specified regions.
#' For how to specify regions, refer to the parameter descriptions and function examples.
#'
#' @param regions_list A vector of regions in the format chrom:start-end.
#' @param regions_bed A bed file with one row for each region you want to determine the CN state from.
#' @param region_names Subset CN states on specific regions (gene symbols e.g FCGR2B).
#' @param these_samples_metadata A metadata table to auto-subset the data to samples in that table before returning.
#' @param all_cytobands Include all cytobands, default is set to FALSE. Currently only supports hg19.
#' @param use_cytoband_name Use cytoband names instead of region names, e.g p36.33.
#'
#' @return Copy number matrix.


#' @import dplyr circlize tibble stringr tidyr
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
                         these_samples_metadata,
                         all_cytobands = FALSE,
                         use_cytoband_name = FALSE){

  this_seq_type="genome" #this only supports genomes currently
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter=this_seq_type)
  }else{
    these_samples_metadata = dplyr::filter(these_samples_metadata,seq_type==this_seq_type)
  }
  if(all_cytobands){
    message("Currently, only grch37 is supported")
  }
  #retrieve the CN value for this region for every segment that overlaps it
  bed2region=function(x){
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }
  if(all_cytobands){
    message("Cytobands are in respect to hg19. This will take awhile but it does work, trust me!")
    use_cytoband_name = TRUE
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
  meta_arranged = these_samples_metadata %>%
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


#' @title GetSample CN Segments.
#'
#' @description Get all segments for a single (or multiple) sample_id(s).
#'
#' @details This function returns CN segments for samples. This works for single sample or multiple samples.
#' For multiple samples, remember to set the Boolean parameter `multiple_samples = TRUE` and give the `sample_lsit` a vector of characters with one sample ID per row.
#' For more information on how this function can be run in different ways, refer to the parameter descriptions, examples and vignettes.
#'
#' @param this_sample_id Optional argument, single sample_id for the sample to retrieve segments for.
#' @param multiple_samples Set to TRUE to return cn segments for multiple samples specified in `samples_list` parameter. Default is FALSE.
#' @param sample_list Optional vector of type character with one sample per row, required if multiple_samples is set to TRUE.
#' @param from_flatfile Set to TRUE by default.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param with_chr_prefix Set to TRUE to add a chr prefix to chromosome names. Default is FALSE.
#' @param streamlined Return a minimal output rather than full details. Default is FALSE.
#'
#' @return A data frame of segments for a specific or multiple sample ID(s).
#'
#' @import dplyr readr RMariaDB DBI
#' @export
#'
#' @examples
#' # Return cn segments for one sample:
#' sample_cn_seg = get_sample_cn_segments(this_sample_id = "some-sample-id", multiple_samples = FALSE)
#'
#' # Return cn segments for multiple samples (provided as vector of sample IDs):
#' samples = get_sample_cn_segments(multiple_samples = TRUE, sample_list = c("some_sample", "another_sample"))
#'
#' # Return cn segments for multiple samples (read csv with one sample per line):
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
    cnv_flatfile_template = check_config_value(config::get("results_flatfiles")$cnv_combined$icgc_dart)
    cnv_path =  glue::glue(cnv_flatfile_template)
    full_cnv_path =  paste0(check_config_value(config::get("project_base")), cnv_path)
    local_full_cnv_path =  paste0(config::get("project_base"), cnv_path)
    if(file.exists(local_full_cnv_path)){
      full_cnv_path = local_full_cnv_path
      #use local file when available
    }
    # check permissions to ICGC data
    permissions = file.access(full_cnv_path, 4)
    if (permissions == -1) {
      message("restricting to non-ICGC data")
      cnv_flatfile_template = check_config_value(config::get("results_flatfiles")$cnv_combined$gambl)
      cnv_path =  glue::glue(cnv_flatfile_template)
      full_cnv_path =  paste0(check_config_value(config::get("project_base")), cnv_path)
    }

    #check for missingness
    if(!file.exists(full_cnv_path)){
      print(paste("missing: ", full_cnv_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    all_segs = suppressMessages(read_tsv(full_cnv_path))
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

      db = check_config_value(config::get("database_name"))
      table_name = check_config_value(config::get("results_tables")$copy_number)
      table_name_unmatched = check_config_value(config::get("results_tables")$copy_number_unmatched)
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

      db = check_config_value(config::get("database_name"))
      table_name = check_config_value(config::get("results_tables")$copy_number)
      table_name_unmatched = check_config_value(config::get("results_tables")$copy_number_unmatched)
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

  #deal with chr prefixes
  if(!with_chr_prefix){
    all_segs = all_segs %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom))
  }else{
    if(!grepl("chr", all_segs$chrom[1])){
      all_segs$chrom = paste0("chr", all_segs$chrom)
      }
  }
  
  if(streamlined){all_segs = dplyr::select(all_segs, ID, CN)}

  return(all_segs)
}


#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range.
#'
#' @details This function returns CN segments for s specified region.
#' There are multiple ways a region can be specified.
#' For example, the user can provide the full region in a "region" format (chr:start-end) to the `region` parameter.
#' Or, the user can provide chromosome, start and end coordinates individually with `chr`, `start`, and `end` parameters.
#' For more usage examples, refer to the parameter descriptions and examples in the vignettes.
#'
#' @param region Region formatted like chrX:1234-5678 or X:1234-56789.
#' @param chromosome The chromosome you are restricting to. Required parameter if region is not specified.
#' @param qstart Start coordinate of the range you are restricting to. Required parameter if region is not specified.
#' @param qend End coordinate of the range you are restricting to. Required parameter if region is not specified.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param this_seq_type Seq type for returned CN segments. Currently, only genome is supported. Capture samples will be added once processed through CNV protocols.
#' @param with_chr_prefix Boolean parameter for toggling if chr prefixes should be present in the return, default is FALSE.
#' @param streamlined Return a basic rather than full MAF format. Default is FALSE.
#' @param from_flatfile Set to TRUE by default.
#'
#' @return A data frame with CN segments for the specified region.
#'
#' @import dplyr readr RMariaDB DBI
#' @export
#'
#' @examples
#' #Example using chromosome, qstart and qend parameters:
#' segments_region_grch37 = get_cn_segments(chromosome = "chr8",
#'                                          qstart = 128723128,
#'                                          qend = 128774067)
#'                                    
#' #Example using the regions parameter:
#' segments_region_hg38 = get_cn_segments(region = "chr8:128,723,128-128,774,067",
#'                                        projection = "hg38",
#'                                        with_chr_prefix = TRUE)
#'
get_cn_segments = function(region,
                           chromosome,
                           qstart,
                           qend,
                           projection = "grch37",
                           this_seq_type = "genome",
                           with_chr_prefix = FALSE,
                           streamlined = FALSE,
                           from_flatfile = TRUE){
  
  #checks
  remote_session = check_remote_configuration(auto_connect = TRUE)
  
  #check seq type and return a message if anything besides "genome" is called. To be updated once capture samples have been processed through CNV protocols.
  if(this_seq_type!="genome"){
    stop("Currently, only genome samples are available for this function. Please select a valid seq type (i.e genome). Compatibility for capture samples will be added soon...")
  }
  
  #get wildcards from this_seq_type (lazy)
  seq_type = this_seq_type
  
  #perform wrangling on the region to have it in the correct format. 
  if(!missing(region)){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }
  
  #deal with chr prefixes for region, based on selected genome projection.
  if(projection == "grch37"){
    if(grepl("chr", chromosome)){
      chromosome = gsub("chr", "", chromosome)
    }
  }else{
    if(!grepl("chr", chromosome)){
      chromosome = paste0("chr", chromosome)
    }
  }
  
  #enforce data type for qend and qstart coordiantes.
  qstart = as.numeric(qstart)
  qend = as.numeric(qend)
  
  if(from_flatfile){
    cnv_flatfile_template = check_config_value(config::get("results_flatfiles")$cnv_combined$icgc_dart)
    cnv_path =  glue::glue(cnv_flatfile_template)
    full_cnv_path =  paste0(check_config_value(config::get("project_base")), cnv_path)
    
    #check permissions to ICGC data.
    permissions = file.access(full_cnv_path, 4)
    if(permissions == -1){
      message("restricting to non-ICGC data")
      cnv_flatfile_template = check_config_value(config::get("results_flatfiles")$cnv_combined$gambl)
      cnv_path =  glue::glue(cnv_flatfile_template)
      full_cnv_path =  paste0(check_config_value(config::get("project_base")), cnv_path)
    }
    
    #check for missingness.
    if(!file.exists(full_cnv_path)){
      print(paste("missing: ", full_cnv_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }
        
    all_segs = suppressMessages(read_tsv(full_cnv_path)) %>%
      dplyr::filter((chrom == chromosome & start <= qstart & end >= qend) | (chrom == chromosome & start >= qstart & end <= qend)) %>%
      as.data.frame()
    
  }else{
    db = check_config_value(config::get("database_name"))
    table_name = check_config_value(config::get("results_tables")$copy_number)
    table_name_unmatched = check_config_value(config::get("results_tables")$copy_number_unmatched)
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

    all_segs_matched = dplyr::tbl(con, table_name) %>%
      dplyr::filter((chrom == chromosome & start <= qstart & end >= qend) | (chrom == chromosome & start >= qstart & end <= qend)) %>%
      as.data.frame() %>%
      dplyr::mutate(method = "battenberg")
      
    all_segs_unmatched = dplyr::tbl(con, table_name_unmatched) %>%
      dplyr::filter((chrom == chromosome & start <= qstart & end >= qend) | (chrom == chromosome & start >= qstart & end <= qend)) %>%
      as.data.frame() %>%
      dplyr::filter(! ID %in% all_segs_matched$ID)  %>%
      dplyr::mutate(method = "controlfreec")
      
      DBI::dbDisconnect(con)
      
      all_segs = rbind(all_segs_matched, all_segs_unmatched)
  }
  
  #mutate CN states.
  all_segs = dplyr::mutate(all_segs, CN = round(2*2^log.ratio))

  #deal with chr prefixes
  if(!with_chr_prefix){
    if(all(str_detect(all_segs$chrom, "chr"))){
      all_segs = all_segs %>%
        dplyr::mutate(chrom = gsub("chr", "", chrom))
    }
  }else{
    if(all(!str_detect(all_segs$chrom, "chr"))){
      all_segs = all_segs %>%
        dplyr::mutate(chrom = paste0("chr", chrom))
    }
  }
  
  #subset to only a few columns with streamlined = TRUE.
  if(streamlined){
    all_segs = dplyr::select(all_segs, ID, CN)
  }
  
  #return data frame with CN segments
  return(all_segs)
}


#' @title Append To Table.
#'
#' @description Housekeeping function to add results to a table.
#'
#' @details INTERNAL FUNCTION, not meant for out-of-package usage.
#'
#' @param table_name The name of the database table to update/populate.
#' @param data_df A dataframe of values to load into the table.
#'
#' @return A table.
#'
#' @import RMariaDB DBI
#'
#' @examples
#' table_up = append_to_table("my_table", "my_df")
#'
append_to_table = function(table_name,
                           data_df){

  db = check_config_value(config::get("database_name"))
  con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  dbWriteTable(con, table_name, table_data, append = TRUE)
}


#' @title Get ASHM Count Matrix.
#'
#' @description Prepare a matrix with one row per sample and one column per region using a set of hypermutated regions.
#'
#' @details Values are the number of mutations in that patient in the region.
#'
#' @param regions_bed A bed file with one row for each region.
#' @param maf_data Optionally provide a data frame in the MAF format, otherwise the database will be used.
#' @param these_samples_metadata This is used to complete your matrix. All GAMBL samples will be used by default. Provide a data frame with at least sample_id for all samples if you are using non-GAMBL data.
#' @param seq_type The seq type to return results for.
#' @param from_indexed_flatfile Boolean parameter set to TRUE per default.
#'
#' @return A matrix.
#'
#' @import dplyr tibble
#' @export
#'
#' @examples
#' regions_bed = grch37_ashm_regions %>% mutate(name = paste(gene, region, sep = "_"))
#' matrix = get_ashm_count_matrix(regions_bed = regions_bed, seq_type="genome")
#'
get_ashm_count_matrix = function(regions_bed,
                                 maf_data,
                                 these_samples_metadata,
                                 seq_type,
                                 from_indexed_flatfile = TRUE){
  if(missing(seq_type)){
    if(missing(these_samples_metadata)){
      stop("Must supply either the seq_type or a metadata data frame from which it can be retrieved")
    }
    seq_type = head(these_samples_metadata) %>% pull(seq_type)
  }
  if(missing(regions_bed)){
    regions_bed = grch37_ashm_regions
  }
  ashm_maf = get_ssm_by_regions(regions_bed = regions_bed,
                                streamlined = TRUE,
                                seq_type=seq_type,
                                maf_data = maf_data,
                                use_name_column = TRUE,
                                from_indexed_flatfile = from_indexed_flatfile)

  ashm_counted = ashm_maf %>%
    group_by(sample_id, region_name) %>%
    tally()

  if(missing(these_samples_metadata)){
    all_meta = get_gambl_metadata(seq_type_filter=seq_type) %>%
      dplyr::select(sample_id)
  }else{
    all_meta = these_samples_metadata %>%
      dplyr::select(sample_id)
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


#' @title Get SSM By Regions.
#'
#' @description Efficiently retrieve all mutations across a range of genomic regions.
#'
#' @details This function internally calls `get_ssm_by_region` to retrieve SSM calls for the specified regions.
#' See parameter descriptions for `get_ssm_by_region` for more information on how the different parameters can be called.
#'
#' @param regions_list Either provide a vector of regions in the chr:start-end format OR.
#' @param regions_bed Better yet, provide a bed file with the coordinates you want to retrieve.
#' @param streamlined Return a basic rather than full MAF format, default is TRUE.
#' @param maf_data Use an already loaded MAF data frame.
#' @param use_name_column If your bed-format data frame has a name column (must be named "name") these can be used to name your regions.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details).
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF
#' @param seq_type The seq_type you want back, default is genome.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param basic_columns Boolean parameter set to FALSE per default. Set to TRUE to return fewer columns.
#'
#' @return Returns a data frame of variants in MAF-like format.
#'
#' @import tibble dplyr tidyr
#' @export
#'
#' @examples
#' #basic usage, adding custom names from bundled ashm data frame
#' regions_bed = grch37_ashm_regions %>% mutate(name = paste(gene, region, sep = "_"))
#' ashm_basic_details = get_ssm_by_regions(regions_bed = regions_bed)
#' full_details_maf = get_ssm_by_regions(regions_bed = regions_bed,basic_columns=T)
#'
get_ssm_by_regions = function(regions_list,
                              regions_bed,
                              streamlined = TRUE,
                              maf_data = maf_data,
                              use_name_column = FALSE,
                              from_indexed_flatfile = TRUE,
                              mode = "slms-3",
                              augmented = TRUE,
                              seq_type = "genome",
                              projection = "grch37",
                              min_read_support = 4,
                              basic_columns = FALSE){


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
    print(regions)
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x,
                                                                streamlined = streamlined,
                                                                from_indexed_flatfile = from_indexed_flatfile,
                                                                mode = mode,
                                                                augmented = augmented,
                                                                seq_type = seq_type,
                                                                projection = projection,
                                                                basic_columns=basic_columns)})
  }else{
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x,
                                                                streamlined = streamlined,
                                                                maf_data = maf_data,
                                                                from_indexed_flatfile = from_indexed_flatfile,
                                                                mode = mode,
                                                                basic_columns=basic_columns)})
  }
  if(!use_name_column){
    rn = regions
  }else{
    rn = regions_bed[["name"]]
  }

  if(basic_columns){
    #this must always force the output to be the standard set.
    #hence, return everything after binding into one data frame
    print("bind_rows")
    return(bind_rows(region_mafs))
  }
  tibbled_data = tibble(region_mafs, region_name = rn)
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)
  if(streamlined){
    unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
      dplyr::select(start, sample_id, region_name)

  }else{
    unlisted_df = mutate(unnested_df, Chromosome = region_mafs$Chromosome, End_Position = region_mafs$End_Position, Start_Position = region_mafs$Start_Position, Tumor_Sample_Barcode = region_mafs$Tumor_Sample_Barcode) %>%
      dplyr::select(Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode, region_name)
  }
  return(unlisted_df)
}


#' @title Get SSM By Region.
#'
#' @description Retrieve all SSMs from the GAMBL database within a single genomic coordinate range.
#'
#' @details This function lets the user specify a region of interest for returning SSM calls within that region.
#' There are multiple ways a region can be specified. For example, the user can provide the full region in a "region" format (chr:start-end) to the `region` parameter.
#' Or, the user can provide chromosome, start and end coordinates individually with `chr`, `start`, and `end` parameters.
#' For more usage examples, refer to the parameter descriptions and examples in the vignettes.
#'
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param qstart Query start coordinate of the range you are restricting to.
#' @param qend Query end coordinate of the range you are restricting to.
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param basic_columns Set to TRUE to override the default behaviour of returning only the first 45 columns of MAF data.
#' @param streamlined Return a basic rather than full MAF format, default is FALSE.
#' @param maf_data Parameter description.
#' @param seq_type The seq_type you want back, default is genome.
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38).
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details).
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF .
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#' @param maf_columns Specify what MAF columns you want back. `basic_columns` needs to be set to TRUE.
#' @param maf_column_types The column types of specified MAF columns `maf_columns`.
#' @param verbose Boolean parameter set to FALSE per default.
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#'
#' @rawNamespace import(vroom, except = c("col_skip", "fwf_positions", "default_locale", "date_names_lang", "cols_only", "output_column", "col_character", "col_guess", "spec", "as.col_spec", "fwf_cols", "cols", "col_date", "col_datetime", "locale", "col_time", "cols_condense", "col_logical", "col_number", "col_integer", "col_factor", "fwf_widths", "date_names_langs", "problems", "date_names", "col_double", "fwf_empty"))
#' @import dplyr RMariaDB DBI stringr
#' @export
#'
#' @examples
#' #basic usage
#' my_mutations = get_ssm_by_region(region = "chr8:128,723,128-128,774,067")
#' #specifying chromosome, start and end individually
#' my_mutations = get_ssm_by_region(chromosome = "8", qstart = 128723128, qend = 128774067)
#' bcl2_all_details = get_ssm_by_region(region="chr18:60796500-60988073",basic_columns=T)
#'
get_ssm_by_region = function(chromosome,
                             qstart,
                             qend,
                             region = "",
                             basic_columns = FALSE,
                             streamlined = FALSE,
                             maf_data,
                             seq_type = "genome",
                             projection = "grch37",
                             from_indexed_flatfile = TRUE,
                             augmented = TRUE,
                             min_read_support = 3,
                             mode = "slms-3",
                             maf_columns = c("Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "t_alt_count"),
                             maf_column_types = "ciici",
                             verbose = FALSE){
  remote_session = check_remote_configuration(auto_connect = TRUE)
  if(basic_columns){
    #this means we ignore/clobber the contents of maf_columns so the first 45 are used instead
    maf_columns = names(maf_header)[c(1:45)]
    maf_column_types = "ccccciiccccccccccccccccccccccnccccccccciiiiii"
    streamlined = FALSE
    #these two arguments are really mutually exclusive so basic_columns must force the other to be FALSE to avoid problems
  }
  #check that maf_columns requested all exist in the header and get their indexes
  if(!all(maf_columns %in% names(maf_header))){
    stop("Cannot find one of the requested maf_columns in your MAF header")
  }
  maf_indexes = maf_header[maf_columns]

  maf_indexes = maf_indexes[order(maf_indexes)]
  maf_columns = names(maf_indexes)
  maf_indexes = unname(maf_indexes)
  #this is to put the indexes and their names back into numerical order because cut returns columns that way

  tabix_bin = check_config_value(config::get("dependencies")$tabix)
  table_name = check_config_value(config::get("results_tables")$ssm)
  db = check_config_value(config::get("database_name"))
  base_path = check_config_value(config::get("project_base"))
  base_path_remote = check_config_value(config::get("project_base",config="default"))

  if(from_indexed_flatfile){

    #test if we have permissions for the full gambl + icgc merge
    if(mode == "slms-3"){
      if(augmented){
        maf_partial_path = check_config_value(config::get("results_flatfiles")$ssm$template$merged$augmented)
      }else{
        maf_partial_path = check_config_value(config::get("results_flatfiles")$ssm$template$merged$deblacklisted)
        }
    }else if (mode == "strelka2"){
      maf_partial_path = check_config_value(config::get("results_flatfiles")$ssm$all$strelka2)
    }else{
      stop("You requested results from indexed flatfile. The mode should be set to either slms-3 (default) or strelka2. Please specify one of these modes.")
    }

    maf_path = glue::glue(maf_partial_path)
    full_maf_path = paste0(base_path, maf_path)
    full_maf_path_comp = paste0(base_path, maf_path, ".bgz")

    if(!file.exists(full_maf_path_comp)){
      print(paste("missing:", full_maf_path_comp))
      check_host(verbose=TRUE)
      if(verbose){
        message("using local file")
        print(paste("HERE:",full_maf_path_comp))
      }
      stop("failed to find the file needed for this")
    }
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
  }else{
    if(projection =="grch37"){
      chromosome = gsub("chr", "", chromosome)
    }
    region=paste0(chromosome,":",qstart,"-",qend)
  }

  if(projection =="grch37"){
    chromosome = gsub("chr", "", chromosome)
  }

 if(missing(maf_data)){
    if(from_indexed_flatfile){
      if(remote_session){
        if(verbose){
          print("ssh session!")
        }
        #Helper function that may come in handy elsewhere so could be moved out of this function if necessary
        run_command_remote = function(ssh_session,to_run){
        output = ssh::ssh_exec_internal(ssh_session,to_run)$stdout
        output = rawToChar(output)
        return(output)
        }

        # NOTE!
        # Retrieving mutations per region over ssh connection is only supporting the basic columns for now in an attempt to keep the transfer of unnecessary data to a minimum

        remote_tabix_bin = check_config_value(config::get("dependencies",config="default")$tabix)

        full_maf_path_comp = paste0(base_path_remote, maf_path, ".bgz")
        #if(!file.exists(full_maf_path_comp)){
        #  message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        #  message('Sys.setenv(R_CONFIG_ACTIVE= "remote")')
        #  check_host()
        #}else{
          message(paste("reading from:", full_maf_path_comp))
        #}

        tabix_command = paste("/home/rmorin/miniconda3/bin/tabix", full_maf_path_comp, region, "| cut -f", paste(maf_indexes,collapse=","))
        if(verbose){
          print(tabix_command)
        }
        #stop()
        muts = run_command_remote(ssh_session,tabix_command)
        muts_region = vroom::vroom(I(muts),col_types = maf_column_types,
                                   col_names=maf_columns)
      }else{

        tabix_command = paste(tabix_bin, full_maf_path_comp, region, "| cut -f" , paste(maf_indexes,collapse=","))
        if(verbose){
          print(tabix_command)
        }
        muts = system(tabix_command, intern = TRUE)
        if(verbose){
          print(paste("TYPES:"))
          print(maf_column_types)
          print("NAMES:")
          print(maf_columns)
        }
        if(length(muts)==0){
          maf_types_sep = str_split(maf_column_types,pattern="")[[1]] %>%
            str_replace_all("c","character") %>%
            str_replace_all("i|n","numeric")

          muts_region = read.table(textConnection(""), col.names = maf_columns,colClasses = maf_types_sep)
        }else{
          muts_region = vroom::vroom(I(muts), col_types = paste(maf_column_types,collapse=""),
                            col_names=maf_columns)
        }
        if(verbose){
          print('SUCCESS')
        }
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
  }

  return(muts_region)
}


#' @title Get Coding SSM.
#'
#' @description Retrieve all coding SSMs from the GAMBL database in MAF-like format.
#'
#' @details Effectively retrieve coding SSM calls. Multiple filtering parameters are available for this function.
#' For more information on how to implement the filtering parameters, refer to the parameter descriptions as well as examples in the vignettes.
#'
#' @param limit_cohort Supply this to restrict mutations to one or more cohorts in a vector.
#' @param exclude_cohort  Supply this to exclude mutations from one or more cohorts in a vector.
#' @param limit_pathology Supply this to restrict mutations to one pathology.
#' @param limit_samples Supply this to restrict mutations to a vector of sample_id (instead of subsetting using the provided metadata)
#' @param these_samples_metadata Supply a metadata table to auto-subset the data to samples in that table before returning
#' @param force_unmatched_samples Optional argument for forcing unmatched samples, using get_ssm_by_samples.
#' @param projection Reference genome build for the coordinates in the MAF file. The default is hg19 genome build.
#' @param seq_type The seq_type you want back, default is genome.
#' @param basic_columns Set to FALSE to override the default behavior of returning only the first 45 columns of MAF data.
#' @param maf_cols if basic_columns is set to FALSE, the user can specify what columns to be returned within the MAF. This parameter can either be a vector of indexes (integer) or a vector of characters (matching columns in MAF).
#' @param from_flatfile Set to TRUE to obtain mutations from a local flatfile instead of the database. This can be more efficient and is currently the only option for users who do not have ICGC data access.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs)
#' @param groups Unix groups for the samples to be included. Default is both gambl and icgc_dart samples.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is TRUE.
#' @param engine Specify one of readr or fread_maf (default) to change how the large files are loaded prior to subsetting. You may have better performance with one or the other but for me fread_maf is faster and uses a lot less RAM.
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#'
#' @import dplyr tidyr RMariaDB DBI
#' @export
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
                          seq_type,
                          basic_columns = TRUE,
                          maf_cols = NULL,
                          from_flatfile = TRUE,
                          augmented = TRUE,
                          min_read_support = 3,
                          groups = c("gambl", "icgc_dart"),
                          include_silent = TRUE,
                          engine = "fread_maf"){

  remote_session = check_remote_configuration()
  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  if(!missing(these_samples_metadata)){
    all_meta = these_samples_metadata
    seq_type = pull(all_meta,seq_type) %>% unique()
    if(length(seq_type)>1){
      stop("More than one seq_type is in this metadata. You can only run this on one seq_type at a time")
    }
  }else{
    if(missing(seq_type)){
      stop("you must provide either seq_type or these_samples_metadata")
    }
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
    maf_template = check_config_value(config::get("results_flatfiles")$ssm$template$cds$deblacklisted)
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(check_config_value(config::get("project_base")), maf_path)
  }

  #get file path for augmented maf
  if(from_flatfile && augmented){
    maf_template = check_config_value(config::get("results_flatfiles")$ssm$template$cds$augmented)
    maf_path = glue::glue(maf_template)
    full_maf_path =  paste0(check_config_value(config::get("project_base")), maf_path)
  }

  #read file
  if(from_flatfile){
    message(paste("reading from:", full_maf_path))

  #check for missingness
    if(!file.exists(full_maf_path)){
      print(paste("missing: ", full_maf_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    if(engine=='fread_maf'){
      if(basic_columns){
        #subset to basic columns during read to save time and memory with lazy loading (in theory)
        select_cols = c(1:45)
        muts = fread_maf(full_maf_path,select_cols=select_cols) %>%
          dplyr::filter(Variant_Classification %in% coding_class) %>%
          as.data.frame()
      }else{
        muts = fread_maf(full_maf_path) %>%
          dplyr::filter(Variant_Classification %in% coding_class) %>%
          as.data.frame()
      }
    }
    mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
    message(paste("mutations from", mutated_samples, "samples"))
  }else{
    #use db if not using flatfile (mostly deprecated)

    table_name = check_config_value(config::get("results_tables")$ssm)
    db = check_config_value(config::get("database_name"))
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



  #subset maf to a specific set of columns (defined in maf_cols)
  if(!is.null(maf_cols) && !basic_columns){
    muts = dplyr::select(muts, all_of(maf_cols))
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


#' @title Get Gene CN and Expression.
#' 
#' @description Get the copy number and expression for a single gene.
#'
#' @details This function works well with both Hugo Symbols and Ensembl Gene IDs. 
#' It's also possible to specify more than one gene.
#'
#' @param gene_symbol One or more gene symbols. Should match the values in a maf file.
#' @param ensembl_id One or more ensembl gene IDs. Only one of hugo_symbols or ensembl_gene_ids may be used.
#'
#' @return A data frame with copy number information and gene expressions.

#' @import dplyr tibble
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


#' @title Get Gene Expression.
#'
#' @description Get the expression for one or more genes for all GAMBL samples.
#'
#' @details Effectively get gene expression for one or multiple genes for al GAMBL samples.
#' This function can also take an already loaded expression matrix (`expression_data`)
#' to prevent the user from having to load the full expression matrix if this function needs to be run in an interactive session.
#' For examples and more info, refer to the parameter descriptions as wella s vignette examples.
#'
#' @param metadata GAMBL metadata.
#' @param hugo_symbols One or more gene symbols. Should match the values in a maf file.
#' @param ensembl_gene_ids One or more ensembl gene IDs. Only one of hugo_symbols or ensembl_gene_ids may be used.
#' @param join_with How to restrict cases for the join. Can be one of genome, mrna or "any".
#' @param all_genes Set to TRUE to return the full expression data frame without any subsetting. Avoid this if you don't want to use tons of RAM.
#' @param expression_data Optional argument to use an already loaded expression data frame (prevent function to re-load full df from flat file or database).
#' @param from_flatfile Deprecated but left here for backwards compatibility.
#'
#' @return A data frame with gene expression.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' MYC_expr = get_gene_expression(hugo_symbols = c("MYC"), join_with = "mrna")
#' # Read full expression values df (no subsetting on genes)
#' full_expression_df = get_gene_expression_new(all_genes = TRUE, join_with = "genome")
#' # Use loaded df (in the previous step) to get expression values for IRF4 and MYC.
#' irf4_myc_expressions = get_gene_expression_new(hugo_symbols = c("IRF4", "MYC"), all_genes = FALSE, join_with = "genome", from_flatfile = FALSE, expression_data = full_expression_df)
#'
get_gene_expression = function(metadata,
                               hugo_symbols,
                               ensembl_gene_ids,
                               join_with = "mrna",
                               all_genes = FALSE,
                               expression_data,
                               from_flatfile = TRUE){

  database_name = check_config_value(config::get("database_name"))
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
  #tidy_expression_file = config::get("results_merged")$tidy_expression_file
  #use combination of base path and relative path instead of full path for flexibility accross sites
  tidy_expression_path = check_config_value(config::get("results_merged")$tidy_expression_path)
  base_path = check_config_value(config::get("project_base"))
  tidy_expression_file = paste0(base_path,tidy_expression_path)
  tidy_expression_file = gsub(".gz$","",tidy_expression_file)

  #check permission and updates paths accordingly
  permissions = file.access(tidy_expression_file, 4)
  if(permissions == -1 ){
    message("restricting to non-ICGC data")
    tidy_expression_path = check_config_value(config::get("results_merged")$tidy_expression_path_gambl)
    tidy_expression_file = paste0(base_path, tidy_expression_path)
  }

  if(!missing(expression_data)){
    tidy_expression_data = as.data.frame(expression_data) #is this necessary? Will it unnecessarily duplicate a large object if it's already a data frame?
    if(!missing(hugo_symbols)){
      #lazily filter on the fly to conserve RAM
      wide_expression_data = tidy_expression_data %>%
        dplyr::filter(Hugo_Symbol %in% hugo_symbols) %>%
        dplyr::select(-ensembl_gene_id) %>%
        group_by(mrna_sample_id,Hugo_Symbol) %>% #deal with non 1:1 mapping of Hugo to Ensembl
        slice_head() %>%
        as.data.frame() %>%
        pivot_wider(names_from = Hugo_Symbol, values_from = expression)
    }else if(!missing(ensembl_gene_ids)){
      wide_expression_data = tidy_expression_data %>%
        dplyr::filter(ensembl_gene_id %in% ensembl_gene_ids) %>%
        dplyr::select(-Hugo_Symbol) %>%
        as.data.frame() %>%
        pivot_wider(names_from = ensembl_gene_id, values_from = expression)
    }else{

      #for when a user wants everything. Need to handle the option of getting back Hugo_Symbol instead
      wide_expression_data = tidy_expression_data %>%
        dplyr::select(-Hugo_Symbol) %>%
        as.data.frame() %>%
        pivot_wider(names_from = ensembl_gene_id, values_from = expression)
    }
  }else{
    if(!file.exists(tidy_expression_file)){
      print(paste("missing: ", tidy_expression_file))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE= "remote")')
      check_host()
    }
    #only ever load the full data frame when absolutely necessary
    if(all_genes & missing(ensembl_gene_ids) & missing(hugo_symbols)){
      wide_expression_data = suppressMessages(read_tsv(tidy_expression_file)) %>%
        as.data.frame() %>%
        pivot_wider(names_from = ensembl_gene_id, values_from = expression)
    }else{
      if(!missing(hugo_symbols)){
        #lazily filter on the fly to conserve RAM (use grep without regex)
        genes_regex=paste(c("-e Hugo_Symbol",hugo_symbols),collapse = " -e ");
        grep_cmd = paste0("grep -w -F ",genes_regex," ",tidy_expression_file)
        print(grep_cmd)
        wide_expression_data = fread(cmd=grep_cmd) %>%
        #wide_expression_data = read_tsv(tidy_expression_file,lazy=TRUE) %>%
          dplyr::select(-ensembl_gene_id) %>%
          dplyr::filter(Hugo_Symbol %in% hugo_symbols) %>%
          group_by(mrna_sample_id,Hugo_Symbol) %>% #deal with non 1:1 mapping of Hugo to Ensembl
          slice_head() %>%
          as.data.frame() %>%
          pivot_wider(names_from = Hugo_Symbol, values_from = expression)
      }
      if(!missing(ensembl_gene_ids)){
        wide_expression_data = suppressMessages(read_tsv(tidy_expression_file,lazy=TRUE)) %>%
          dplyr::select(-Hugo_Symbol) %>%
          dplyr::filter(ensembl_gene_id %in% ensembl_gene_ids) %>%
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


#' @title Get Manta SV By Samples.
#'
#' @description Load the manta output for a set of samples.
#'
#' @details This is a convenience wrapper function for get_manta_sv_by_sample (and called by get_manta_sv).
#'
#' @param these_samples_metadata The only required parameter is a metadata table (data frame) that must contain a row for each sample you want the data from. The additional columns the data frame needs to contain, besides sample_id, are: unix_group, genome_build, seq_type, pairing_status.
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Default value is 0.1.
#' @param min_score The lowest Manta somatic score for a SV to be returned. Default value is 40.
#' @param pass If set to TRUE, only return SVs that are annotated with PASS in the FILTER column. Set to FALSE to keep all variants, regardless if they PASS the filters. Default is TRUE. 
#' @param projection The projection of returned calls. Default is grch37.
#' 
#' @return A data frame containing the Manta outputs from all sample_id in these_samples_metadata in a bedpe-like format with additional columns extracted from the VCF column.
#'
#' @import dplyr stringr
#' @export
#'
#' @examples
#' all_sv = get_manta_sv
#' missing_samples = get_gambl_metadata() %>% anti_join(all_sv, by = c("sample_id" = "tumour_sample_id"))
#' missing_from_merge = get_manta_sv_by_samples(these_samples_metadata = missing_samples)
#'
get_manta_sv_by_samples = function(these_samples_metadata,
                                   min_vaf = 0.1,
                                   min_score = 40,
                                   pass = TRUE,
                                   projection = "grch37"){
  
  #check remote configuration
  remote_session = check_remote_configuration(auto_connect = TRUE)

  #get sample IDs from metadata.
  samples = pull(these_samples_metadata, sample_id)
  
  #create an empty list.
  all_bedpe = list()
  
  #wrap get_manta_sv_by_sample.
  all_bedpe = lapply(samples, function(x){get_manta_sv_by_sample(this_sample_id = x,
                                                                 these_samples_metadata = these_samples_metadata,
                                                                 force_lift = FALSE, #the wrapper function performs liftover on all samples that need it.
                                                                 return_anyway = TRUE, #make sure unlifted calls, with the extra column (need_lift) are returned.
                                                                 min_vaf = min_vaf,
                                                                 min_score = min_score,
                                                                 pass = pass,
                                                                 projection = projection)})
  
  #un-nest list into long format.
  merged_bedpe = bind_rows(all_bedpe)
  
  #take out all calls that need to be lifted
  to_be_lifted = merged_bedpe %>%
    dplyr::filter(need_lift == TRUE)
  
  #lift to selected projection
  lifted_calls = liftover_bedpe(bedpe_df = to_be_lifted, target_build = projection)
  
  #subset calls that does not need a "lift"
  no_lift_needed = merged_bedpe %>%
    dplyr::filter(need_lift == FALSE)
  
  #combine calls (lifted and not lifted), arrange and sort accordingly, drop temporary column
  merged_bedpe = bind_rows(lifted_calls, no_lift_needed) %>%
    dplyr::select(-need_lift)

  #add chr prefix to the chromosome name for builds that expect it, but only add when necessary
  #hg38 and hg19 (with chr prefix)
  if(projection %in% c("hg38", "hg19")){
    merged_bedpe = merged_bedpe %>%
      dplyr::mutate(CHROM_A = case_when(str_detect(CHROM_A, "chr") ~ CHROM_A, TRUE ~ paste0("chr", CHROM_A))) %>%
      dplyr::mutate(CHROM_B = case_when(str_detect(CHROM_B, "chr") ~ CHROM_B, TRUE ~ paste0("chr", CHROM_B)))
  }

  #grch37 and grch38 (no chr prefix)
  if(projection %in% c("grch37", "grch38")){
    merged_bedpe = merged_bedpe %>%
      dplyr::mutate(CHROM_A = gsub("chr", "", CHROM_A)) %>%
      dplyr::mutate(CHROM_B = gsub("chr", "", CHROM_B))
  }

  #sort data frame
  merged_bedpe = merged_bedpe %>%
    arrange(CHROM_A, CHROM_B, START_A)

  #return merged manta SVs.
  return(merged_bedpe)
}


#' @title Get Manta SV By Sample.
#'
#' @description Load the manta output (from individual flat file) for 1 sample.
#'
#' @details This function is used for retrieving Manta results (structural variants) from individual flat-files (one sample). 
#' For multiple samples, please see `get_manta_sv_by_samples` (a convenience wrapper function for `get_manta_by_sample`). 
#' Additional columns are extracted from the VCF column and standard filtering options are available. 
#' This function also performs a lift-over to selected projection, if needed. 
#' Please note, if `force_lift` is set to FALSE, an extra column will be added that states if the returned variant calls need to be lifted. 
#' The value for this column is returned TRUE (for all rows) if the available genome projection for the selected sample does not match the selected projection (i.e requiring the user to manually lift the calls).
#'
#' @param this_sample_id The single sample ID you want to obtain the result from. If this parameter is not supplied, the function will retrieve sample ID from the supplied metadata table (these_samples_metadata).
#' @param these_samples_metadata A metadata table containing metadata for this_sample_id, or sample of interest. This parameter is required.
#' @param force_lift If TRUE, coordinates will be lifted (if needed) to the selected projection. Default is FALSE. WARNING: if your code calls this function directly, set this parameter to TRUE to ensure that the returned calls are in respect to the requested projection.
#' @param return_anyway Set to TRUE to force variant calls to be returned, even if they're not lifted, This parameter should only ever be modified from the default setting when this function is called by another function that handles the liftOver separately.
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Default value is 0.1.
#' @param min_score The lowest Manta somatic score for a SV to be returned. Default value is 40.
#' @param pass If set to TRUE, only return SVs that are annotated with PASS in the FILTER column. Set to FALSE to keep all variants, regardless if they PASS the filters. Default is TRUE. 
#' @param projection The projection of returned calls. Default is grch37.
#'
#' @return a data frame containing the Manta outputs from this_sample_id in a bedpe-like format with additional columns extracted from the VCF column.
#'
#' @import config dplyr readr stringr tibble
#' @export
#'
#' @examples
#' #example 1
#' #get manta calls for a sample that needs to be lifted to "hg38" and let this function take care of the liftover step for you. 
#' my_sv = get_manta_sv_by_sample(this_sample_id = "99-27783_tumorA", these_samples_metadata = get_gambl_metadata(), projection = "hg38", force_lift = TRUE)
#'
#' #example 2
#' #get manta calls based on an already filtered metadata (with one sample ID)
#' my_metadata = get_gambl_metadata() %>% dplyr::filter(sample_id=="99-27783_tumorA")
#' my_sv = get_manta_sv_by_sample(these_samples_metadata = my_metadata, projection = "hg38", force_lift = TRUE)
#' 
get_manta_sv_by_sample = function(this_sample_id,
                                  these_samples_metadata,
                                  force_lift = FALSE,
                                  return_anyway = FALSE,
                                  min_vaf = 0.1,
                                  min_score = 40,
                                  pass = TRUE,
                                  projection = "grch37"){

  #safetynet for preventing users to mistakenly return un-lifted variant calls.
  if(!force_lift){ #i.e I will run liftover on my own, based on the information in the extra column (need_lift).
    if(!return_anyway){ 
      stop("If you know what you are doing and wish to liftover the returned sample yourself, set return_anyway to TRUE. If you want this function to handle the liftover for you, set force_lift = TRUE")
    }
  }
  
  #check remote configuration
  remote_session = check_remote_configuration(auto_connect = TRUE)
  
  if(missing(this_sample_id)){
    if(!nrow(these_samples_metadata) == 1){
      stop("There is more than one sample in the supplied metadata table. Either subset metadata to only have one sample, provide the this_sample_id parameter OR consider running get_manta_sv_by_samples")
    }
    this_sample_id = these_samples_metadata$sample_id
  }
  
  these_samples_metadata = dplyr::filter(these_samples_metadata, sample_id == this_sample_id)
  if(!nrow(these_samples_metadata==1)){
    stop("metadata does not seem to contain your this_sample_id or you didn't provide one")
  }
  
  #get wildcards
  tumour_sample_id = this_sample_id
  unix_group = pull(these_samples_metadata, unix_group)
  seq_type = pull(these_samples_metadata, seq_type)
  genome_build = pull(these_samples_metadata, genome_build)
  pairing_status = pull(these_samples_metadata, pairing_status)
  
  if(pairing_status == "matched"){
    normal_sample_id = pull(these_samples_metadata, normal_sample_id)
  }else{
    normal_sample_id = config::get("unmatched_normal_ids")[[unix_group]][[seq_type]][[genome_build]]
  }

  #get samples from individual flat files
  path_template = check_config_value(config::get("results_flatfiles")$sv_manta$template)
  
  if(!remote_session){
    path_template_full = paste0(check_config_value(config::get("project_base")), path_template)
    bedpe_path = glue::glue(path_template_full)
    if(!file.exists(bedpe_path)){
      print(paste("missing: ", bedpe_path))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }
  }else{
    local_path_template = paste0(check_config_value(config::get("project_base", config = "remote")), path_template)
    bedpe_path = glue::glue(local_path_template)

    #check if the requested file is on your local machine, if not, get it!
    if(!file.exists(bedpe_path)){
      remote_path_template = paste0(check_config_value(config::get("project_base", config = "default")), path_template)
      remote_bedpe_path = glue::glue(remote_path_template)
      cat(paste0("Local file not found.\ntrying to copy requested file: ", remote_bedpe_path, "\n", "To: ", bedpe_path))
      dirN = dirname(bedpe_path)
      suppressMessages(suppressWarnings(dir.create(dirN, recursive = T)))
      ssh::scp_download(ssh_session, remote_bedpe_path, dirN)
    }
  }
  
  #read sample flat-file
  message(paste0("Reading ", this_sample_id, " from: ", bedpe_path))
  bedpe_dat_raw = suppressMessages(read_tsv(bedpe_path, comment = "##", col_types = "cddcddccccccccccccccccc"))

  #return empty data frame
  if(!nrow(bedpe_dat_raw==0)){
    message(paste0("WARNING! No SV calls found in flat-file for: ", this_sample_id))
    return()
  }
  
  #if the selected projection is different from the genome build (for the selected sample), add information that this sample needs to be lifted (by get_manta_by_samples).
  if(genome_build != projection){
    if(all(str_detect(genome_build, "37|19"))){
      if(projection %in% c("grch37", "hg19")){
        bedpe_dat_raw = bedpe_dat_raw %>%
          add_column(need_lift = FALSE)
      }else if(projection %in% c("hg38", "grch38")){
        bedpe_dat_raw = bedpe_dat_raw %>%
          add_column(need_lift = TRUE)
      }else{
        stop(paste0(projection, " is not a valid projection, acceptable projections are; grch37, grch38, hg19, hg38"))
      }
    }else if(all(str_detect(genome_build, "38"))){
      if(projection %in% c("grch37", "hg19")){
        bedpe_dat_raw = bedpe_dat_raw %>%
          add_column(need_lift = TRUE)
      }else if(projection %in% c("hg38", "grch38")){
        bedpe_dat_raw = bedpe_dat_raw %>%
          add_column(need_lift = FALSE)
      }else{
        stop(paste0(projection, " is not a valid projection, acceptable projections are; grch37, grch38, hg19, hg38"))
      }
    }
  }else{
    bedpe_dat_raw = bedpe_dat_raw %>%
      add_column(need_lift = FALSE)
  }
  
  if(force_lift){
    if(bedpe_dat_raw$need_lift[1] == TRUE){
      bedpe_dat_raw = liftover_bedpe(bedpe_df = bedpe_dat_raw, target_build = projection)
      message(paste0(this_sample_id, " flat-file is not available in the selected projection, running liftover_bedpe..."))
      message(paste0(this_sample_id, " successfully lifted to ", projection)) 
    }
  }

  #data wrangling
  #get infos
  infos = pull(bedpe_dat_raw, tumour_sample_id)
  infos_n = pull(bedpe_dat_raw, normal_sample_id)

  #create new columns with sample IDs
  bedpe_dat = bedpe_dat_raw %>%  
    mutate(tumour_sample_id = tumour_sample_id, normal_sample_id = normal_sample_id, pair_status = pairing_status)

  #rename columns to match the expected format
  colnames(bedpe_dat)[c(1:6)] = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B")

  #extract info fields from VCF
  #tumour sample
  bedpe_dat$VAF_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
  bedpe_dat$DP_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")),2)[1])})

  #normal sample
  bedpe_dat$VAF_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
  bedpe_dat$DP_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 2)[1])})

  #get somatic score
  bedpe_dat$SCORE = sapply(bedpe_dat$INFO_A, function(x){as.numeric(tail(unlist(strsplit(x, "=")), 1))})

  #Rename and select columns to match what is returned with get_combined_sv.
  bedpe_dat = bedpe_dat %>%
    rename("DP" = "DP_tumour", "manta_name" = "ID") %>%
    dplyr::select("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B",
                  "manta_name", "SCORE", "STRAND_A", "STRAND_B", "tumour_sample_id",
                  "normal_sample_id", "VAF_tumour", "DP", "pair_status", "FILTER", "need_lift")
  
  #VAF and somatic score filtering.
  bedpe_dat = bedpe_dat %>%
    dplyr::filter(VAF_tumour >= min_vaf & SCORE >= min_score)
  
  #Filter on FILTER (variant callers variant filter criteria).
  if(pass){
    bedpe_dat = bedpe_dat %>%
      dplyr::filter(FILTER == "PASS")
  }
  
  #Deal with chr prefixes based on projection
  if(force_lift){
    #hg38 and hg19 (with chr prefix)
    if(projection %in% c("hg38", "hg19")){
      bedpe_dat = bedpe_dat %>%
        dplyr::mutate(CHROM_A = case_when(str_detect(CHROM_A, "chr") ~ CHROM_A, TRUE ~ paste0("chr", CHROM_A))) %>%
        dplyr::mutate(CHROM_B = case_when(str_detect(CHROM_B, "chr") ~ CHROM_B, TRUE ~ paste0("chr", CHROM_B)))
    }
    
    #grch37 and grch38 (no chr prefix)
    if(projection %in% c("grch37", "grch38")){
      bedpe_dat = bedpe_dat %>%
        dplyr::mutate(CHROM_A = gsub("chr", "", CHROM_A)) %>%
        dplyr::mutate(CHROM_B = gsub("chr", "", CHROM_B))
    }
    
    #remove the additional column (need_lift)
    bedpe_dat = bedpe_dat %>%
      dplyr::select(-need_lift)
  }
  
  #enforce column types and sort returned calls
  bedpe_dat = bedpe_dat %>%
    mutate(across(c(CHROM_A, CHROM_B, manta_name, STRAND_A, STRAND_B, tumour_sample_id, normal_sample_id, pair_status, FILTER), as.character)) %>%
    mutate(across(c(START_A, END_A, START_B, END_B, SCORE, VAF_tumour, DP), as.numeric)) %>%
    arrange(CHROM_A, CHROM_B, START_A)

  return(bedpe_dat)
}
