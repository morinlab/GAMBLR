#require packages
require("dbplyr")
require("tidyverse")
require("data.table")

#' Initialize a new cBioPortal instance or update existing portal data set, can also be used to retrieve sample ids included in study.
#'
#' @param short_name A concise name for your portal project.
#' @param include_icgc_data Whether or not you want ICGC and other external data included.
#' @param human_friendly_name A slightly more verbose name for your project.
#' @param project_name Unique ID for your project.
#' @param description A verbose description of your data set.
#' @param cancer_type Cancer types included in study, default is "mixed".
#' @param data_cancer_type The cancer type abbreviation, e.g., "brca". (for data_cancer_type file).
#' @param data_cancer_name  The name of the cancer type, e.g., "Breast Invasive Carcinoma". (for data_cancer_type file).
#' @param data_dedicated_colour CSS color name of the color associated with this cancer study, e.g., "HotPink". (for data_cancer_type file).
#' @param data_cancer_parent_type The type_of_cancer field of the cancer type of which this is a subtype, e.g., "Breast". You can set parent to tissue, which is the reserved word to place the given cancer type at "root" level in the "studies oncotree" that will be generated in the homepage (aka query page) of the portal. (for data_cancer_type file).
#' @param exclude_cohorts Cohorts to be excluded.
#' @param gambl_maf maf origin.
#' @param gambl_icgc_maf Icgc maf origin.
#' @param overwrite Flag to specify that files should be overwritten if they exist. Default is TRUE.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return A vector of sample_id for the patients that have been included.
#' @export
#'
#' @examples
#'
#' Setup study and save included ids as a vector of characters:
#' ids = setup_study(out_dir = "GAMBLR/cBioPortal/instance01/")
#'
setup_study = function(short_name = "GAMBL",
                       include_icgc_data = FALSE,
                       human_friendly_name = "GAMBL data",
                       project_name = "gambl_minus_icgc",
                       description = "GAMBL data without ICGC",
                       cancer_type = "mixed",
                       data_cancer_type = "brca", 
                       data_cancer_name = "Invasive Breast Carcinoma", 
                       data_dedicated_colour = "HotPink", 
                       data_cancer_parent_type = "Breast",
                       exclude_cohorts = c("FFPE_Benchmarking"),
                       gambl_maf = "maf_slms3_hg19",
                       gambl_icgc_maf = "maf_slms3_hg19_icgc",
                       overwrite = TRUE,
                       out_dir){
  
  #determine what table to query and what restrictions to use for the MAF data
  if(include_icgc_data){
    maf_table = gambl_icgc_maf
  }else{
    maf_table = gambl_maf
  }
  
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
  
  #meta cancer type
  meta_cancer_type = paste0(out_dir, "meta_cancer_type.txt")
  
  meta_cancer_type_content = paste0("genetic_alteration_type: CANCER_TYPE\n", 
                                    "datatype: CANCER_TYPE\n", 
                                    "data_filename: data_cancer_type.txt\n")
  
  cat(meta_cancer_type_content, file = meta_cancer_type)
  
  #data cancer type
  data_cancer_type = paste0(out_dir, "data_cancer_type.txt")
  
  data_cancer_type_content = paste(data_cancer_type, 
                                   data_cancer_name, 
                                   data_dedicated_colour, 
                                   data_cancer_parent_type, "\n",
                                   sep = "\t")
  
  cat(data_cancer_type_content, file = data_cancer_type)
  
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
                                  "profile_description: Mutation data from whole genome sequencing.\n", 
                                  "profile_name: Mutations\n", 
                                  "data_filename: ", data_mutations, "\n", 
                                  "swissprot_identifier: name\n")
  
    cat(meta_mutations_content, file = meta_mutations)
  
  if(overwrite){
    #create the actual MAF file by querying the database using the API
    coding_ssms = get_coding_ssm(exclude_cohort = exclude_cohorts)
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


#' Annotate SVs and create the input for fusions to be displayed in cBioPortal instance.
#'
#' @param short_name A concise name for your portal project.
#' @param include_icgc_data Whether or not you want ICGC and other external data included.
#' @param human_friendly_name A slightly more verbose name for your project.
#' @param project_name Unique ID for your project.
#' @param gambl_maf maf origin.
#' @param gambl_icgc_maf Icgc maf origin.
#' @param description A verbose description of your data set.
#' 
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return A vector of sample_id for the patients that have been included.
#' @export
#'
#' @examples
#'
#' Setup fusions and save included sample ids into a vector of character
#' fusion_ids = setup_fusions(out_dir = "GAMBLR/cBioPortal/instance01/)
#' 
setup_fusions = function(short_name = "GAMBL",
                         include_icgc_data = FALSE,
                         human_friendly_name = "GAMBL data",
                         project_name = "gambl_minus_icgc",
                         description = "GAMBL data without ICGC",
                         gambl_maf = "maf_slms3_hg19",
                         gambl_icgc_maf = "maf_slms3_hg19_icgc",
                         out_dir){

  #determine what table to query and what restrictions to use for the MAF data
  if(include_icgc_data){
    maf_table = gambl_icgc_maf
  }else{
    maf_table = gambl_maf
  }

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
  unannotated_sv = get_manta_sv()

  annotated_sv = annotate_sv(unannotated_sv) %>%
    dplyr::filter(!is.na(partner)) %>%
    as.data.frame()

  fusion_samples = pull(annotated_sv, tumour_sample_id) %>%
    unique()

  #deal with any cases not in metadata
  fusions_df =  data.frame(Hugo_Symbol = annotated_sv$gene,
                           Entrez_Gene_Id = annotated_sv$entrez,
                           Center = "BCGSC",
                           Tumor_Sample_Barcode = annotated_sv$tumour_sample_id,
                           Fusion = c(pull(unite(annotated_sv, fusion, partner, gene, sep = "-"), fusion)),
                           DNA_support = "yes",
                           RNA_support = "no",
                           Method = "Manta",
                           Frame = "in-frame")

  fusions_df = distinct(fusions_df, Tumor_Sample_Barcode, Fusion, .keep_all = TRUE)

  #determine what table to query and what restrictions to use for the MAF data
  if(include_icgc_data){
    maf_table = gambl_icgc_maf
  }else{
    maf_table = gambl_maf
  }

  nfkbiz_entrez = 64332
  
  nfkbiz_utr_ssm = get_ssm_by_gene(gene_symbol = "NFKBIZ") %>%
    dplyr::filter(Variant_Classification == "3'UTR") %>%
    pull(Tumor_Sample_Barcode) %>%
    unique()

  nfkbiz.mut.df = data.frame(Hugo_Symbol = "NFKBIZ",
                             Entrez_Gene_Id = nfkbiz_entrez,
                             Center = "BCGSC",
                             Tumor_Sample_Barcode = nfkbiz_utr_ssm,
                             Fusion = "NFKBIZ-UTR",
                             DNA_support = "yes",
                             RNA_support = "no",
                             Method = "SLMS-3",
                             Frame = "in-frame")
  
  #get any SV breakpoints that are in the 3'UTR of NFKBIZ
  nfkbiz_utr_region = "chr3:101,578,185-101,579,902"
  data_fusions = paste0(out_dir, "data_fusions.txt")

  nfkbiz.svs = get_manta_sv(region = nfkbiz_utr_region) %>%
    pull(tumour_sample_id) %>%
    unique()

  nfkbiz.sv.df = data.frame(Hugo_Symbol = "NFKBIZ",
                            Entrez_Gene_Id = nfkbiz_entrez,
                            Center = "BCGSC",
                            Tumor_Sample_Barcode = nfkbiz.svs,
                            Fusion = "NFKBIZ-SV",
                            DNA_support = "yes",
                            RNA_support = "no",
                            Method = "Manta",
                            Frame = "in-frame")

  all_fusions = rbind(fusions_df, nfkbiz.sv.df, nfkbiz.mut.df)
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


#' Finish setting up a new cBioPortal instance or updating an existing portal data set.
#'
#' @param short_name A concise name for your portal project.
#' @param include_icgc_data Whether or not you want ICGC and other external data included.
#' @param human_friendly_name A slightly more verbose name for your project.
#' @param project_name Unique ID for your project.
#' @param description A verbose description of your data set.
#' @param cancer_type Cancer types included in study, default is "mixed".
#' @param sample_ids A vector of all the sample_id that were included in any of the data files for cbioportal.
#' @param overwrite Flag to specify that files should be overwritten if they exist. Default is TRUE.
#' @param out_dir The full path to the base directory where the files are being created.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' Finalize study using generated list of ids from setup_study and setup_fusions:
#' finalize_study(sample_ids = c(ids, fusion_ids), out_dir = "GAMBLR/cBioPortal/instance01/)
#'
finalize_study = function(short_name = "GAMBL",
                          include_icgc_data = FALSE,
                          human_friendly_name = "GAMBL data",
                          project_name = "gambl_minus_icgc",
                          description = "GAMBL data without ICGC",
                          cancer_type = "mixed",
                          sample_ids,
                          overwrite = TRUE,
                          out_dir){
  
  #create necessary files
  #create case list
  caselist = paste0(out_dir, "case_lists/cases_sequenced.txt")
  
  tabseplist = paste(unique(sample_ids), collapse = "\t")
  
  caselistdata = c(paste0("cancer_study_identifier: ", project_name), 
                   paste0("stable_id: ", project_name, "_sequenced"), "case_list_name: Samples sequenced.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.", 
                   paste0(c("case_list_ids:", tabseplist),collapse = " "))
  
  cat(caselistdata, sep = "\n", file = caselist)

  #create case list all
  caselist_all = paste0(out_dir, "case_lists/cases_all.txt")
  
  caselistdata = c(paste0("cancer_study_identifier: ", project_name), 
                   paste0("stable_id: ",project_name, "_allcases"), "case_list_name: Samples sequenced.", "case_list_description: This is this case list that contains all samples that are profiled for mutations.", 
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))
  
  cat(caselistdata, sep = "\n", file = caselist_all)

  #meta samples
  #prepare and write out the relevant metadata
  clinsamp = paste0(out_dir, "data_clinical_samples.txt")
  meta_samples = get_gambl_metadata() %>%
    dplyr::filter(sample_id %in% sample_ids) %>%
    dplyr::select(patient_id, sample_id, pathology, EBV_status_inf, cohort, time_point, ffpe_or_frozen, myc_ba, bcl6_ba, bcl2_ba, COO_consensus, DHITsig_consensus, lymphgen)

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


#' Helper function for checking integrity of study files.
#'
#' @param data_clinical_samples Path to clinical file.
#' @param data_fusions Path to data_fusion file from setup_fusions.
#' @param cases_fusions Path to cases_fusion from setup_fusions.
#' @param cases_all Path to cases_all from setup_study.
#' @param cases_sequenced Path to cases_sequenced from setup_study.
#' @param project_name Project name, should match what is specified under setup_study/setup_fusions.
#' @param out_dir Directory with all study related files, the only argument that needs to be specified, given that path to all generated study files are not changed from default.
#' 
#' @return Nothing.
#' @export
#'
#' @examples
#' Perform sanity check on generated files (i.e ensure that sample ids in case lists are represented in clinical file and filter out such samples, if any):
#' samples_not_in_clinical = study_check(out_dir = "GAMBLR/cBioPortal/instance01/)
#'
study_check = function(data_clinical_samples_path = "data_clinical_samples.txt", 
                       data_fusions_path = "data_fusions.txt",
                       cases_fusions_path = "case_lists/cases_fusion.txt",
                       cases_all_path = "case_lists/cases_all.txt",
                       cases_sequenced_path = "case_lists/cases_sequenced.txt",
                       project_name = "gambl_minus_icgc",
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
