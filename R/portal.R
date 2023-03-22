#' @title Setup Study (cBioPortal).
#'
#' @description Initialize a new cBioPortal instance or update existing portal data set, can also be used to retrieve sample ids included in study.
#'
#' @details This function internally calls `get_coding_ssm` to retrieve coding mutations to be included in the study (if `overwrite = TRUE`).
#' In addition, this function also creates and sets up the proper folder hierarchy and writes the files necessary to import a new cBioPortal study.
#' Before a study is ready to be imported to cBioPortal, the user also needs to run `setup_fusions` and `finalize_study`.
#' Optionally the user can also run `study_check` to ensure all samples described by the "clinical" file are included in the study.
#' Also, note that the parameters chosen for this function have to match the same parameters called for any subsequent study function calls.
#'
#' @param seq_type_filter the seq type you are setting up a study for, default is "genome".
#' @param short_name A concise name for your portal project.
#' @param human_friendly_name A slightly more verbose name for your project.
#' @param project_name Unique ID for your project.
#' @param description A verbose description of your data set.
#' @param overwrite Flag to specify that files should be overwritten if they exist. Default is TRUE.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return A vector of sample_id for the patients that have been included.
#'
#' @import dplyr readr
#' @export
#'
#' @examples
#' Setup study and save included ids as a vector of characters:
#' \dontrun{
#' ids = setup_study(out_dir = "GAMBLR/cBioPortal/instance01/")
#' }
#' 
setup_study = function(seq_type_filter = "genome",
                       short_name = "GAMBL",
                       human_friendly_name = "GAMBL data",
                       project_name = "gambl_genome",
                       description = "GAMBL data from genome",
                       overwrite = TRUE,
                       out_dir){
  cancer_type="mixed"

  #set up the new directory
  if (!file.exists(out_dir)){
    dir.create(out_dir)
    dir.create(paste0(out_dir, "case_lists"))
  }else{
    dir.create(paste0(out_dir, "case_lists"))
  }

  #create necessary files
  #meta study
  meta_study = paste0(out_dir, "meta_study.txt")

  meta_study_content = paste0("type_of_cancer: ", cancer_type, "\n",
                              "cancer_study_identifier: ", project_name, "\n",
                              "name: ", human_friendly_name, "\n",
                              "short_name: ", short_name, "\n",
                              "description: ", description, "\n",
                              "add_global_case_list: true\n")

  cat(meta_study_content, file = meta_study)


  #meta clinical samples
  meta_clinical_samples = paste0(out_dir, "meta_clinical_samples.txt")

  meta_clinical_samples_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                         "genetic_alteration_type: CLINICAL\n",
                                         "datatype: SAMPLE_ATTRIBUTES\n",
                                         "data_filename: data_clinical_samples.txt")

  cat(meta_clinical_samples_content, file = meta_clinical_samples)

  #meta clinical patients
  meta_clinical_patients = paste0(out_dir, "meta_clinical_patient.txt")

  meta_clinical_patients_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                          "genetic_alteration_type: CLINICAL\n",
                                          "datatype: PATIENT_ATTRIBUTES\n",
                                          "data_filename: data_clinical_patient.txt")

  cat(meta_clinical_patients_content, file = meta_clinical_patients)

  #meta mutations
  meta_mutations = paste0(out_dir, "meta_mutations_extended.txt")
  data_mutations = "data_mutations_extended.maf"

  meta_mutations_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                  "genetic_alteration_type: MUTATION_EXTENDED\n",
                                  "datatype: MAF\n",
                                  "stable_id: mutations\n",
                                  "show_profile_in_analysis_tab: true\n",
                                  "profile_description: Mutation data\n",
                                  "profile_name: Mutations\n",
                                  "data_filename: ", data_mutations, "\n",
                                  "swissprot_identifier: name\n")

  cat(meta_mutations_content, file = meta_mutations)

  if(overwrite){
    #create the actual MAF file by querying the database using the API
    coding_ssms = get_coding_ssm(seq_type = seq_type_filter)
    data_mutations_full = paste0(out_dir, "data_mutations_extended.maf")
    write_tsv(coding_ssms, data_mutations_full, na = "")
  }else{
    #read in the MAF instead
    coding_ssms = data.table::fread(file = data_mutations_full,
                                    sep = "\t",
                                    stringsAsFactors = FALSE,
                                    verbose = FALSE,
                                    data.table = TRUE,
                                    showProgress = TRUE,
                                    header = TRUE,
                                    fill = TRUE,
                                    skip = "Hugo_Symbol",
                                    quote = "")
  }
  ids = coding_ssms %>%
    pull(Tumor_Sample_Barcode) %>%
    unique()

  return(ids)
}


#' @title Setup Fusions (cBioPortal).
#'
#' @description Annotate SVs and create the input for fusions to be displayed in cBioPortal instance.
#'
#' @details This function calls `get_combined_sv` and runs `annotate_sv` on the returned data frame.
#' Should be run as the next step after running `setup_study`. Note that the parameters called with this function
#' has to match the parameter arguments of `setup_study`, i.e if `short_name` is for `setup_study` is "GAMBL",
#' then the `short_name` in `setup_fusions` also has to be "GAMBL", etc.
#'
#' @param short_name A concise name for your portal project.
#' @param human_friendly_name A slightly more verbose name for your project.
#' @param project_name Unique ID for your project.
#' @param gambl_maf maf origin.
#' @param gambl_icgc_maf ICGC maf origin.
#' @param description A verbose description of your data set.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return A vector of sample_id for the patients that have been included.
#'
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' fusion_ids = setup_fusions(out_dir = "GAMBLR/cBioPortal/instance01/")
#' }
setup_fusions = function(short_name = "GAMBL",
                         human_friendly_name = "GAMBL data",
                         project_name = "gambl_genome",
                         description = "GAMBL data from genome",
                         gambl_maf = "maf_slms3_hg19",
                         gambl_icgc_maf = "maf_slms3_hg19_icgc",
                         out_dir){


  #create necessary files
  #meta fusions
  meta_fusions = paste0(out_dir, "meta_fusions.txt")

  meta_fusion_content = paste0("cancer_study_identifier: ", project_name, "\n",
                               "genetic_alteration_type: FUSION\n",
                               "datatype: FUSION\n",
                               "stable_id: fusion\n",
                               "show_profile_in_analysis_tab: true\n",
                               "profile_name: Fusions\n",
                               "profile_description: Fusion data\n",
                               "data_filename: data_fusions.txt\n")

  cat(meta_fusion_content, file = meta_fusions)

  #get SV breakpoints and annotate them
  unannotated_sv = get_combined_sv()

  annotated_sv = annotate_sv(unannotated_sv) %>%
    dplyr::filter(!is.na(partner)) %>%
    as.data.frame()

  fusion_samples = pull(annotated_sv, tumour_sample_id) %>%
    unique()

  #deal with any cases not in metadata
  fusions_df =  data.frame(Hugo_Symbol = annotated_sv$gene,
                           Center = "BCGSC",
                           Tumor_Sample_Barcode = annotated_sv$tumour_sample_id,
                           Fusion = c(pull(unite(annotated_sv, fusion, partner, gene, sep = "-"), fusion)),
                           DNA_support = "yes",
                           RNA_support = "no",
                           Method = "SVAR",
                           Frame = "in-frame")

  fusions_df = distinct(fusions_df, Tumor_Sample_Barcode, Fusion, .keep_all = TRUE)

  #determine what table to query and what restrictions to use for the MAF data

  # TO DO: Fix this code to work with the indexed MAF file using get_ssm_by_region instead of by gene
  #nfkbiz_entrez = 64332
  #nfkbiz_utr_ssm = get_ssm_by_gene(gene_symbol = "NFKBIZ") %>%
  #  dplyr::filter(Variant_Classification == "3'UTR") %>%
  #  pull(Tumor_Sample_Barcode) %>%
  #  unique()

  #nfkbiz.mut.df = data.frame(Hugo_Symbol = "NFKBIZ",
  #                           Entrez_Gene_Id = nfkbiz_entrez,
  #                           Center = "BCGSC",
  #                           Tumor_Sample_Barcode = nfkbiz_utr_ssm,
  #                           Fusion = "NFKBIZ-UTR",
  #                           DNA_support = "yes",
  #                           RNA_support = "no",
  #                           Method = "SLMS-3",
  #                           Frame = "in-frame")

  #get any SV breakpoints that are in the 3'UTR of NFKBIZ
  #nfkbiz_utr_region = "chr3:101,578,185-101,579,902"
  data_fusions = paste0(out_dir, "data_fusions.txt")
  #TODO: FIX! this is also broken
  #nfkbiz.svs = get_combined_sv(region = nfkbiz_utr_region) %>%
  #  pull(tumour_sample_id) %>%
  #  unique()

  #nfkbiz.sv.df = data.frame(Hugo_Symbol = "NFKBIZ",
  #                          Entrez_Gene_Id = nfkbiz_entrez,
  #                          Center = "BCGSC",
  #                          Tumor_Sample_Barcode = nfkbiz.svs,
  #                          Fusion = "NFKBIZ-SV",
  #                          DNA_support = "yes",
  #                          RNA_support = "no",
  #                          Method = "Manta",
  #                          Frame = "in-frame")

  #all_fusions = rbind(fusions_df, nfkbiz.sv.df, nfkbiz.mut.df)
  all_fusions = fusions_df
  fusion.cases = as.character(unique(all_fusions$Tumor_Sample_Barcode))
  write_tsv(all_fusions, data_fusions)

  #create necessary files
  #create caselist fusions
  caselist_fusion = paste0(out_dir, "case_lists/cases_fusion.txt")

  tabseplist = paste(fusion.cases, collapse = "\t")

  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, "_fusions"), "case_list_name: Samples with fusions.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))

  cat(caselistdata, sep = "\n", file = caselist_fusion)
  return(fusion.cases)
}


#' @title Finalize Study (cBioPortal).
#'
#' @description Finish setting up a new cBioPortal instance or updating an existing portal data set.
#'
#' @details This function should be run as the last (or third step) in setting up a new cBioPortal instance.
#' The functions that should be run prior to these functions are; `setup_study` and `setup_fusions`.
#' `finalize_study` creates all the necessary tables and metadata files (case lists) that are required to import a new study into cBioPortal.
#' Note, that all parameter arguments used in this function have to match the same parameter arguments for the previously run functions (`setup_study` and `setup_fusions`).
#'
#' @param seq_type_filter the seq type you are setting up a study for, default is "genome".
#' @param short_name A concise name for your portal project.
#' @param human_friendly_name A slightly more verbose name for your project.
#' @param project_name Unique ID for your project.
#' @param description A verbose description of your data set.
#' @param cancer_type Cancer types included in study, default is "mixed".
#' @param these_sample_ids A vector of all the sample_id that were included in any of the data files for cBioPortal (i.e the output from `setup_study` and `setup_fusions`).
#' @param overwrite Flag to specify that files should be overwritten if they exist. Default is TRUE.
#' @param meta_columns Optional parameter for specifying metadata fields in the study. If not provided, the function will resort to using default metadata columns.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return Nothing.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' finalize_study(these_sample_ids = c(ids, fusion_ids), out_dir = "GAMBLR/cBioPortal/instance01/")
#' }
finalize_study = function(seq_type_filter = "genome",
                          short_name = "GAMBL",
                          human_friendly_name = "GAMBL data",
                          project_name = "gambl_genome",
                          description = "GAMBL data from genome",
                          cancer_type = "mixed",
                          these_sample_ids,
                          overwrite = TRUE,
                          meta_columns,
                          out_dir){

  #resort to "default" metadata columns, if not specified
  if(missing(meta_columns)){
    meta_columns = c("patient_id", "sample_id", "pathology",
                     "EBV_status_inf", "cohort", "time_point",
                     "ffpe_or_frozen", "myc_ba", "bcl6_ba",
                     "bcl2_ba", "COO_consensus", "DHITsig_consensus", "lymphgen")
  }

  #create necessary files
  #create case list
  caselist = paste0(out_dir, "case_lists/cases_sequenced.txt")

  tabseplist = paste(unique(these_sample_ids), collapse = "\t")

  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, "_sequenced"), "case_list_name: Samples sequenced", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                   paste0(c("case_list_ids:", tabseplist),collapse = " "))

  cat(caselistdata, sep = "\n", file = caselist)

  #create case list all
  caselist_all = paste0(out_dir, "case_lists/cases_all.txt")

  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ",project_name, "_allcases"), "case_list_name: Samples sequenced", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))

  cat(caselistdata, sep = "\n", file = caselist_all)

  #meta samples
  #prepare and write out the relevant metadata
  clinsamp = paste0(out_dir, "data_clinical_samples.txt")
  meta_samples = get_gambl_metadata(seq_type_filter=seq_type_filter) %>%
    dplyr::filter(sample_id %in% these_sample_ids) %>%
    dplyr::select(all_of(meta_columns))

  colnames(meta_samples) = toupper(colnames(meta_samples))

  header = paste0("#Patient Identifier\tSample Identifier\tSubtype\tEBV status\tCohort\tTime point\tFFPE\tMYC_BA\tBCL6_BA\tBCL2_BA\tCOO\tDHITsig\tLymphGen\n",
                  "#Patient identifier\tSample Identifier\tSubtype\tEBV status\tCohort\tTime point\tFFPE\tMYC_BA\tBCL6_BA\tBCL2_BA\tCOO\tDHITsig\tLymphGen\n",
                  "#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\n",
                  "#1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\n")

  cat(header, file = clinsamp)

  write.table(meta_samples, file = clinsamp, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)

  #create clinical meta data
  #first get the patient_id list
  clinpat = paste0(out_dir, "data_clinical_patient.txt")
  patient_ids = pull(meta_samples, PATIENT_ID)

  all_outcomes = get_gambl_outcomes(time_unit = "month", censor_cbioportal = TRUE, patient_ids = patient_ids, complete_missing = TRUE) %>%
    dplyr::select(c("patient_id", "OS_STATUS", "OS_MONTHS", "DFS_STATUS", "DFS_MONTHS", "age", "sex"))

  colnames(all_outcomes) = toupper(colnames(all_outcomes))

  header = paste0("#Patient Identifier\tOverall Survival Status\tOverall Survival (Months)\tDisease Free Status\tDisease Free (Months)\tAGE\tSEX\n",
                  "#Patient Identifier\tOverall Survival Status\tOverall Survival (Months)\tDisease Free Status\tDisease Free (Months)\tAge\tSex\n",
                  "#STRING\tSTRING\tNUMBER\tSTRING\tNUMBER\tNUMBER\tSTRING\n",
                  "#1\t1\t1\t1\t1\t1\t1\n")

  cat(header, file = clinpat)

  write.table(all_outcomes, file = clinpat, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
}


#' @title Study Check (cBioPortal).
#'
#' @description Helper function for checking integrity of study files.
#'
#' @details This function was designed to ensure that all the sample IDs described in the maf are actually present in the clinical files.
#' If this is not the case, the function will notify the user what samples are found in the case list that are not described in the clinical file.
#' The function then sub-sets the case list to only include samples from the clinical file.
#' Note that the `project_name` has to match what is specified for the previously run functions (i.e `setup_study`, `setup_fusions` and `finalize_study`).
#'
#' @param data_clinical_samples_path Path to clinical file.
#' @param data_fusions_path Path to data_fusion file from setup_fusions.
#' @param cases_fusions_path Path to cases_fusion from setup_fusions.
#' @param cases_all_path Path to cases_all from setup_study.
#' @param cases_sequenced_path Path to cases_sequenced from setup_study.
#' @param project_name Project name, should match what is specified under setup_study/setup_fusions.
#' @param out_dir Directory with all study related files, the only argument that needs to be specified, given that paths to all generated study files are not changed from default.
#'
#' @return Nothing.
#'
#' @import dplyr readr
#' @export
#'
#' @examples
#' \dontrun{
#' samples_not_in_clinical = study_check(out_dir = "GAMBLR/cBioPortal/instance01/")
#' }
study_check = function(data_clinical_samples_path = "data_clinical_samples.txt",
                       data_fusions_path = "data_fusions.txt",
                       cases_fusions_path = "case_lists/cases_fusion.txt",
                       cases_all_path = "case_lists/cases_all.txt",
                       cases_sequenced_path = "case_lists/cases_sequenced.txt",
                       project_name = "gambl_genome",
                       out_dir){

  #read clinical file (skip header)
  data_clinical_samples = data.table::fread(file = paste0(out_dir, data_clinical_samples_path), sep = "\t", header = FALSE, skip = 5)

  #cases fusions
  cases_fusion = data.table::fread(file = paste0(out_dir, cases_fusions_path), sep = "	", skip = 4, header = FALSE)

  #transform data and strip irrelevant characters
  cases_fusion = t(cases_fusion) %>%
    as.data.frame() %>%
    mutate_at("V1", str_replace, "case_list_ids: ", "")

  #return samples that are present in case_lists but not in clinical file
  not_in_case_lists = setdiff(cases_fusion$V1, data_clinical_samples$V2)

  if(length(not_in_case_lists) > 0){
    #print message
    message("Warning! Samples found in case lists that are not described in the clinical file (data_clincaal_samples.txt). The following samples will be removed from all case lists and data fusions:")
    print(not_in_case_lists)

    #intersect with clinical file to remove sample ids that are not represented (in clinical file)
    cases_fusions_check = intersect(cases_fusion$V1, data_clinical_samples$V2)

    #put updated sample ids back into case list and overwrite original file
    caselist_fusions = paste0(out_dir, cases_fusions_path)

    tabseplist_fusions = paste(cases_fusions_check, collapse = "\t")

    caselistdata_fusions = c(paste0("cancer_study_identifier: ", project_name),
                             paste0("stable_id: ", project_name, "_sequenced"), "case_list_name: Samples sequenced.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                             paste0(c("case_list_ids:", tabseplist_fusions),collapse = " "))

    cat(caselistdata_fusions, sep = "\n", file = paste0(out_dir, cases_fusions_path))

    #data fusions
    data_fusions = data.table::fread(file = paste0(out_dir, data_fusions_path), sep = "\t", header = TRUE)

    #remove all samples not in clinical file that are present in data fusions
    data_fusions = data_fusions[!grepl(paste(not_in_case_lists, collapse = "|"), data_fusions$Tumor_Sample_Barcode),]

    #overwrite data_fusions.txt
    data_fusions_out = paste0(out_dir, data_fusions_path)
    write_tsv(data_fusions, data_fusions_out)

    #cases all
    cases_all = data.table::fread(file = paste0(out_dir, cases_sequenced_path), sep = "	", skip = 4, header = FALSE)

    #transform data and strip irrelevant characters
    cases_all = t(cases_all) %>%
      as.data.frame() %>%
      mutate_at("V1", str_replace, "case_list_ids: ", "")

    #intersect with clinical file to remove sample ids that are not represented (in clinical file)
    cases_all_check = intersect(cases_all$V1, data_clinical_samples$V2)

    #put updated sample ids back into case list and overwrite original file
    caselist_all = paste0(out_dir, cases_all_path)

    tabseplist_all = paste(cases_all_check, collapse = "\t")

    caselistdata_all = c(paste0("cancer_study_identifier: ", project_name),
                         paste0("stable_id: ", project_name, "_sequenced"), "case_list_name: Samples sequenced.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                         paste0(c("case_list_ids:", tabseplist_all),collapse = " "))

    cat(caselistdata_all, sep = "\n", file = paste0(out_dir, cases_all_path))

    #cases sequenced
    cases_sequenced = data.table::fread(file = paste0(out_dir, cases_sequenced_path), sep = "	", skip = 4, header = FALSE)

    #transform data and strip irrelevant characters
    cases_sequenced = t(cases_sequenced) %>%
      as.data.frame() %>%
      mutate_at("V1", str_replace, "case_list_ids: ", "")

    #intersect with clinical file to remove sample ids that are not represented (in clinical file)
    cases_sequenced_check = intersect(cases_sequenced$V1, data_clinical_samples$V2)

    #put updated sample ids back into case list and overwrite original file
    caselist_sequenced = paste0(out_dir, cases_sequenced_path)

    tabseplist_sequenced = paste(cases_sequenced_check, collapse = "\t")

    caselistdata_sequenced = c(paste0("cancer_study_identifier: ", project_name),
                               paste0("stable_id: ", project_name, "_sequenced"), "case_list_name: Samples sequenced.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                               paste0(c("case_list_ids:", tabseplist_sequenced),collapse = " "))

    cat(caselistdata_sequenced, sep = "\n", file = paste0(out_dir, cases_sequenced_path))
  }else{
    message("No additional samples were found in case lists that are not described in the clinical file (data_clinical_samples.txt)")
  }
}


#' @title Custom cBioPortal case list.
#'
#' @description Create and a custom case list for easy data subset in cBioPortal.
#'
#' @details Convenience function for specifying custom case lists that can be browsed on cBioPortal.
#' This function takes a set of sample IDs `these_sample_ids` and intersect the IDs with what's available in the study-specific clinical file.
#' This function also extracts the project name for the specified study, i.e the project name that is defined withing the folder specified under the `dir` parameter.
#'
#' @param these_sample_ids A vector of sample IDs to be subset into a case list. Required parameter.
#' @param caselist_name Name of the generated case list (name does not include the file format, this will be added automagically). This parameter is required.
#' @param caselist_description A verbose description of the created case list. Required.
#' @param return_missing_samples Boolean parameter. Set to TRUE to return all sample IDs that are in the desired case list, but not represented in the study specific clinical file. Default is FALSE.
#' @param dir The directory where all study specific files live.
#'
#' @return
#' 
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' #get some sample IDs
#' my_samples = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL") %>%
#'  dplyr::filter(cohort == "FL_GenomeCanada") %>%
#'  pull(sample_id)
#'
#' #create case list with selected sample IDs
#' custom_caselist(these_sample_ids = my_samples,
#'                 caselist_name = "FL_Canada",
#'                 caselist_description = "Follicular Lymphoma from the Genome Canada Study",
#'                 dir = "../path/to/study_directory/")
#' }
#' 
custom_caselist = function(these_sample_ids,
                           caselist_name,
                           caselist_description,
                           return_missing_samples = FALSE,
                           dir){
  
  #get path to the clinical file holding all sample IDs
  clinical_file = data.table::fread(file = paste0(dir, "data_clinical_samples.txt"), sep = "\t", header = FALSE, skip = 5)
  
  #get project name
  meta_study = data.table::fread(file = paste0(dir, "meta_study.txt"), fill = TRUE)
  project_name = pull(meta_study[2,2])
  
  #pull sample IDs
  clinical_ids = clinical_file %>%
    pull(V2)
  
  #intersect sample IDs from the clinical file with specified sample IDs (these_sample_ids) to ensure no IDs are described in the new case list, that are not represented in the clinical file.
  ids = intersect(clinical_ids, these_sample_ids)
  
  #print how many sample, if any, were removed in the process.
  not_in_clin = setdiff(these_sample_ids, clinical_ids)
  if(!return_missing_samples){
    print(paste0(length(not_in_clin), " Samples not described in the clinical file (data.clinical_samples.txt) for the selected study. To see what samples are not represented, set return_missing_samples = TRUE."))
  }else{
    print(paste0(length(not_in_clin), " Samples not described in the clinical file (data.clinical_samples.txt) for the selected study. The missing samples are: "))
    print(not_in_clin)
  }
  
  #create new case list file
  new_caselist = paste0(dir, "case_lists/", caselist_name, ".txt")
  
  #get sample ID subset belonging to the defined case list
  tabseplist = paste(unique(ids), collapse = "\t")
  
  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, caselist_name),
                   paste0("case_list_name: ", caselist_name),
                   paste0("case_list_description: ", caselist_description),
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))
  
  cat(caselistdata, sep = "\n", file = new_caselist)
}


#' @title Get Study Info.
#'
#' @description Function for retrieving study specific identifiers.
#'
#' @details This function takes one required parameter (`dir`).
#' This is the relative path to the main directory for the study of interest.
#' The function reads in the meta_study.txt file and extracts unique study identifiers.
#' This returns a list that holds all the identifiers.
#' The user can also return the study identifiers to the global environment.
#' To do so, set the `list_to_global = TRUE`.
#'
#' @param dir The relative path to the study directory (expects to find meta_study.txt in this folder).
#' @param list_to_global Boolean parameter, if set to TRUE all study identifiers will be returned to the global environment. Default is FALSE.
#'
#' @return A list with study specific identifiers, or nothing (i.e list_to_global = TRUE).
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' #return study identifiers as  list:
#' my_study_info= get_study_info(dir = "path/to/study/")
#'
#' #return all identifiers to the global environment:
#' get_study_info(dir = "path/to/study/", list_to_global = TRUE)
#' }
#' 
get_study_info = function(dir,
                          list_to_global = FALSE){

  #check if meta data file exists in the selected directory
  if(file.exists(paste0(dir, "meta_study.txt"))){
    study_meta =  data.table::fread(file = paste0(dir, "meta_study.txt"), sep = "\t", header = FALSE)
  }else{
    stop("Unable to find meta_study.txt in the specified folder (dir)...")
  }

  #tidy the data frame
  meta_info = gsub(".*: ","", study_meta$V1)

  #create a list with study identifiers
  meta_list = list(cancer_type = meta_info[1],
                   project_name = meta_info[2],
                   human_friendly_name = meta_info[3],
                   short_name = meta_info[4],
                   description = meta_info[5])

  #add vectors to global environment
  if(list_to_global){
    list2env(meta_list, envir = .GlobalEnv)
  }else{
    return(meta_list)
  }
}


#' @title Create cBioPortal Study.
#'
#' @description Wrapper function for creating a import-ready cBioPortal study.
#'
#' @details This function internally calls `setup_study`, `setup_fusions`, `finalize_study` and `study_check` to generate all necessary files for importing a study into cBioPortal.
#' This function was developed to streamline this step and at the same time ensure that the study information and selected data type is consistent throughout the individual steps of generating a study.
#' In addition, the user can also control if the generated study should be checked for sample IDs in case lists that are not described in the clinical file.
#' This potentially will prevent an annoying error that prevents the study to be imported into the active cBioPortal instance, default is TRUE.
#' Fusions are also handled based on the selected seq type (`this_seqtype`).
#'
#' @param this_seqtype The seq type you want to generate a study for. Default is "genome".
#' @param short_name A concise name for your portal project. Default is "GAMBL".
#' @param human_friendly_name A slightly more verbose name for your project. Default is "GAMBL data".
#' @param project_name Unique ID for your project. Default is "gambl_all".
#' @param description A verbose description of your data set. This is what the study will be named when accessing it through cBioPortal. Default is "GAMBL data from genome".
#' @param gambl_maf MAF origin.
#' @param gambl_icgc_maf ICGC MAF origin.
#' @param cancer_type Cancer types included in study, default is "mixed".
#' @param overwrite Flag to specify that files should be overwritten if they exist. Default is TRUE.
#' @param check_study Boolean parameter that controls if the generated study should be checked for sample IDs in case lists, that are not described in the clinical file. Default is TRUE.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return Nothing. Rather, this function generates all files necessary for successfully importing a study into an active cBioPortal instance.
#'
#' @examples
#' \dontrun{
#' #generate cBioPortal study for all GAMBL genome samples:
#' cbioportal_create()
#'
#' #generate a cBioPortal study for all GAMVL capture samples:
#' cbioportal_create(this_seqtype = "capture", description = "GAMBL data from exomes")
#' }
#' 
cbioportal_create = function(this_seqtype = "genome",
                             short_name = "GAMBL",
                             human_friendly_name = "GAMBL data",
                             project_name = "gambl_genome",
                             description = "GAMBL data from genome",
                             gambl_maf = "maf_slms3_hg19",
                             gambl_icgc_maf = "maf_slms3_hg19_icgc",
                             cancer_type = "mixed",
                             overwrite = TRUE,
                             check_study = TRUE,
                             out_dir){

  #setup study
  ids = setup_study(seq_type_filter = this_seqtype,
                    short_name = short_name,
                    human_friendly_name = human_friendly_name,
                    project_name = project_name,
                    description = description,
                    overwrite = overwrite,
                    out_dir = out_dir)

  #setup fusions
  if(this_seqtype == "genome"){
    fusion_ids = setup_fusions(short_name = short_name,
                               human_friendly_name = human_friendly_name,
                               project_name = project_name,
                               description = description,
                               gambl_maf = gambl_maf,
                               gambl_icgc_maf = gambl_icgc_maf,
                               out_dir = out_dir)

  }else if(this_seqtype == "capture"){
    message("The selected seq type is capture, no fusions will be generated...")

  }else{
    stop("Please enter a valid seq_type (i.e genome or capture)...")
  }

  #finalize study
  if(this_seqtype == "genome"){
    finalize_study(seq_type_filter = this_seqtype,
                   short_name = short_name,
                   human_friendly_name = human_friendly_name,
                   project_name = project_name,
                   description = description,
                   cancer_type = cancer_type,
                   these_sample_ids = c(ids, fusion_ids),
                   overwrite = overwrite, out_dir = out_dir)

  }else if(this_seqtype == "capture"){
    finalize_study(seq_type_filter = this_seqtype,
                   short_name = short_name,
                   human_friendly_name = human_friendly_name,
                   project_name = project_name,
                   description = description,
                   cancer_type = cancer_type,
                   these_sample_ids = ids,
                   overwrite = overwrite, out_dir = out_dir)

  }else{
    stop("Please enter a valid seq_type (i.e genome or capture)...")
  }

  #check study
  study_check(data_clinical_samples_path = "data_clinical_samples.txt",
              data_fusions_path = "data_fusions.txt",
              cases_fusions_path = "case_lists/cases_fusion.txt",
              cases_all_path = "case_lists/cases_all.txt",
              cases_sequenced_path = "case_lists/cases_sequenced.txt",
              project_name = project_name,
              out_dir = out_dir)
}


#' @title Setup Expression Data (cBioPortal).
#'
#' @description Generate expression data based on a set of genes, format and export data for cBioPortal.
#'
#' @details This function takes a set of genes with the `these_genes` (character of vectors) parameter and returns expression data.
#' Expression data is then formatted to match the expected format for import to a cBioPortal study.
#' If no genes are provided, the function will default to all genes that are defined in the `lymphoma_genes` bundled data.
#' This function internally calls `get_gene_expression` for returning expression data as outlined above.
#'
#' @param project_name Unique ID for your project.
#' @param clinical_file_path The path to the study specific clinical file (data_clinical_samples.txt).
#' @param these_genes Specify a set of genes (character of vectors) that you want to return expression data for. If no genes are provided, this function will resort to all lymphoma genes.
#' @param expression_df Optional argument for providing an already loaded expression matrix.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return A vector of characters with sample IDs that expression data was generated for.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' #return expression data for lymphoma genes (all samples)
#' expression_ids = setup_expression_data(out_dir = "../../")
#' }
setup_expreession_data = function(project_name = "gambl_genome",
                                  clinical_file_path = "data_clinical_samples.txt",
                                  these_genes,
                                  expression_df,
                                  out_dir){
  
  #get path to the clinical file holding all sample IDs
  clinical_file = data.table::fread(file = paste0(out_dir, clinical_file_path), sep = "\t", header = FALSE, skip = 5)
  
  #pull sample IDs
  clinical_ids = clinical_file %>%
    pull(V2)
  
  #default to all lymphoma genes if no specific genes are supplied for getting expression values.
  if(missing(these_genes)){
    these_genes = lymphoma_genes %>% pull(Gene)
  }
  
  if(missing(expression_df)){
    #get expression data
    expression_matrix = get_gene_expression(hugo_symbols = these_genes, join_with = "genome")
  }else{
    expression_matrix = expression_df
  }
  
  #expression metadata
  meta_expression = paste0(out_dir, "meta_expression.txt")
  
  meta_expression_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                   "genetic_alteration_type: MRNA_EXPRESSION\n",
                                   "datatype: CONTINUOUS\n",
                                   "stable_id: rna_seq_mrna\n",
                                   "show_profile_in_analysis_tab: false\n",
                                   "profile_name: mRNA expression (microarray)\n",
                                   "profile_description: Expression levels (Alignment microarray)\n",
                                   "data_filename: data_expressions.txt\n")
  
  cat(meta_expression_content, file = meta_expression)
  
  #subset sample IDs from the expression matrix
  expression_samples = pull(expression_matrix, sample_id) %>%
    unique()
  
  #intersect sample IDs from the clinical file with specified sample IDs (these_sample_ids) to 
  #ensure no IDs are described in the new case list, that are not represented in the clinical file.
  ids = intersect(clinical_ids, expression_samples)
  
  #filter expression matrix on the sample IDs that are defined in the clinical file
  expression_matrix = expression_matrix %>% dplyr::filter(sample_id %in% ids) %>%
    dplyr::select(-capture_sample_id)
  
  #transform the expression matrix to match expected format
  tmp_exp <- data.frame(t(expression_matrix[]))
  names(tmp_exp) <- tmp_exp[1,]
  tmp_exp <- tmp_exp[-1,]
  expression_matrix <- tibble::rownames_to_column(tmp_exp, "Hugo_Symbol")
  
  #write expression matrix to file
  data_expression_full = paste0(out_dir, "data_expressions.txt")
  write_tsv(df3, data_expression_full)
  
  #create case list for expression data
  caselist_expression = paste0(out_dir, "case_lists/cases_expression.txt")
  
  tabseplist = paste(ids, collapse = "\t")
  
  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, "_expression"),
                   paste0("case_list_name: ", "Expression Data"),
                   paste0("case_list_description: ", "This is this case list that contains all samples that have expression data."),
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))
  
  cat(caselistdata, sep = "\n", file = caselist_expression)
  
  return(ids)
}
