require("dbplyr")
require("tidyverse")
require("data.table")

#functions for creating a cBioportal instance using GAMBL data
#some global variables that we will probably change later
gambl_db = "gambl_test"
gambl_maf = "maf_slms3_hg19"
gambl_icgc_maf = "maf_slms3_hg19_icgc"


#' Annotate SVs and create the input for fusions to be displayed in cBioPortal instance
#'
#' @param short_name a concise name for your portal project
#' @param include_icgc_data whether or not you want ICGC and other external data included
#' @param human_friendly_name a slightly more verbose name for your project
#' @param project_name unique ID for your project
#' @param description A verbose description of your data set
#' @param out_dir The full path to the base directory where the files are being created
#'
#' @return a vector of sample_id for the patients that have been included
#' @export
#'
#' @examples
setup_fusions = function(short_name = "GAMBL",
                         include_icgc_data = FALSE,
                         human_friendly_name = "GAMBL data",
                         project_name = "gambl_minus_icgc",
                         description = "GAMBL data without ICGC",
                         out_dir){
  meta_fusions = paste0(out_dir, "meta_fusions.txt")
  data_fusions = paste0(out_dir, "data_fusions.txt")
  fusions_detailed = paste0(out_dir, "annotated_fusions_detail.tsv")
  caselist_fusion = paste0(out_dir, "case_lists/cases_fusion.txt")

  #determine what table to query and what restrictions to use for the MAF data
  #TODO: fix this once we have the ICGC SV data in the database
  if(include_icgc_data){
    maf_table = gambl_icgc_maf
  }else{
    maf_table = gambl_maf
  }

  meta_fusion_content = paste0("cancer_study_identifier: ", project_name,"\n",
                               "genetic_alteration_type: FUSION\n",
                               "datatype: FUSION\n",
                               "stable_id: fusion\n",
                               "show_profile_in_analysis_tab: true\n",
                               "profile_name: Fusions\n",
                               "profile_description: Fusion data\n",
                               "data_filename: data_fusions.txt\n"

  )
  cat(meta_fusion_content, file = meta_fusions)

  #get SV breakpoints and annotate them

  unannotated_sv = get_manta_sv() #no filters

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
  nfkbiz_utr_ssm = get_ssm_by_gene(table = maf_table, gene_symbol = "NFKBIZ") %>%
    dplyr::filter(Variant_Classification == "3'UTR") %>% 
    pull(Tumor_Sample_Barcode) %>% 
    unique()

  nfkbiz.mut.df = data.frame(Hugo_Symbol = "NFKBIZ",
                             Entrez_Gene_Id = nfkbiz_entrez,
                             Center = "BCGSC",
                             Tumor_Sample_Barcode = nfkbiz_utr_ssm,
                             Fusion = "NFKBIZ-UTR",
                             DNA_support = "yes",
                             RNA_support="no",
                             Method = "SLMS-3",
                             Frame = "in-frame")
  #get any SV breakpoints that are in the 3'UTR of NFKBIZ
  nfkbiz_utr_region = "chr3:101,578,185-101,579,902"

  nfkbiz.svs= get_manta_sv(region=nfkbiz_utr_region) %>% pull(tumour_sample_id) %>% unique()

  nfkbiz.sv.df = data.frame(Hugo_Symbol = "NFKBIZ",
                            Entrez_Gene_Id = nfkbiz_entrez,
                            Center = "BCGSC",
                            Tumor_Sample_Barcode = nfkbiz.svs,
                            Fusion = "NFKBIZ-SV",
                            DNA_support = "yes",
                            RNA_support="no",
                            Method = "Manta",
                            Frame = "in-frame")

  all_fusions = rbind(fusions_df, nfkbiz.sv.df, nfkbiz.mut.df)

  fusion.cases = as.character(unique(all_fusions$Tumor_Sample_Barcode))

  write_tsv(all_fusions, data_fusions)

  tabseplist = paste(fusion.cases, collapse = "\t")
  caselistdata = c(paste0("cancer_study_identifier: ",project_name),
                   paste0("stable_id: ", project_name, "_fusions"),
                   "case_list_name: Samples with fusions.",
                   "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))

  cat(caselistdata, sep = "\n", file = caselist_fusion)

  return(fusion.cases)
}


#' Finish setting up a new cBioPortal instance or updating an existing portal data set
#'
#' @param short_name a concise name for your portal project
#' @param include_icgc_data whether or not you want ICGC and other external data included
#' @param human_friendly_name a slightly more verbose name for your project
#' @param project_name unique ID for your project
#' @param description A verbose description of your data set
#' @param out_dir The full path to the base directory where the files are being created
#' @param sample_ids A vector of all the sample_id that were included in any of the data files for cbioportal
#'
#' @return nothing
#' @export
#'
#' @examples
finalize_study = function(short_name = "GAMBL",
                          include_icgc_data = FALSE,
                          human_friendly_name = "GAMBL data",
                          project_name = "gambl_minus_icgc",
                          description = "GAMBL data without ICGC",
                          cancer_type = "mixed",
                          out_dir,
                          overwrite = FALSE,
                          sample_ids){
  caselist = paste0(out_dir, "case_lists/cases_sequenced.txt")
  caselist_all = paste0(out_dir, "case_lists/cases_all.txt")
  clinsamp = paste0(out_dir, "data_clinical_samples.txt")
  clinpat = paste0(out_dir, "data_clinical_patient.txt")

  tabseplist = paste(unique(sample_ids), collapse = "\t")
  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, "_sequenced"),
                   "case_list_name: Samples sequenced.",
                   "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))
  cat(caselistdata, sep = "\n", file = caselist)

  caselistdata = c(paste0("cancer_study_identifier: ", project_name),
                   paste0("stable_id: ", project_name, "_allcases"),
                   "case_list_name: Samples sequenced.",
                   "case_list_description: This is this case list that contains all samples that are profiled for mutations.",
                   paste0(c("case_list_ids:", tabseplist), collapse = " "))
  cat(caselistdata, sep = "\n", file = caselist_all)

  #prepare and write out the relevant metadata

  meta_samples = get_gambl_metadata() %>%
    dplyr::filter(sample_id %in% sample_ids) %>%
    dplyr::select(patient_id, sample_id, pathology, EBV_status_inf, cohort, time_point, ffpe_or_frozen, myc_ba, bcl6_ba, bcl2_ba, COO_consensus, DHITsig_consensus, lymphgen)

  colnames(meta_samples) = toupper(colnames(meta_samples))

  header=paste0("#Patient Identifier\tSample Identifier\tSubtype\tEBV status\tCohort\tTime point\tFFPE\tMYC_BA\tBCL6_BA\tBCL2_BA\tCOO\tDHITsig\tLymphGen\n",
                "#Patient identifier\tSample Identifier\tSubtype\tEBV status\tCohort\tTime point\tFFPE\tMYC_BA\tBCL6_BA\tBCL2_BA\tCOO\tDHITsig\tLymphGen\n",
                "#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\n",
                "#1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\n")
  cat(header,file=clinsamp)
  write.table(meta_samples,file=clinsamp,sep="\t",row.names=F,quote=F,append = TRUE)

  #get the clinical metadata
  #first get the patient_id list
  patient_ids = pull(meta_samples,PATIENT_ID)

  all_outcomes = get_gambl_outcomes(time_unit="month",censor_cbioportal = TRUE,patient_ids=patient_ids,complete_missing=TRUE) %>%
    dplyr::select(c("patient_id","OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS","age","sex"))
  colnames(all_outcomes) = toupper(colnames(all_outcomes))

  header=paste0("#Patient Identifier\tOverall Survival Status\tOverall Survival (Months)\tDisease Free Status\tDisease Free (Months)\tAGE\tSEX\n",
                "#Patient Identifier\tOverall Survival Status\tOverall Survival (Months)\tDisease Free Status\tDisease Free (Months)\tAge\tSex\n",
                "#STRING\tSTRING\tNUMBER\tSTRING\tNUMBER\tNUMBER\tSTRING\n",
                "#1\t1\t1\t1\t1\t1\t1\n")

  cat(header, file = clinpat)
  write.table(all_outcomes, file = clinpat, sep = "\t", row.names = F, quote = F, append = TRUE)

}


#' Initialize a new cBioPortal instance or update existing portal data set
#'
#' @param short_name a concise name for your portal project
#' @param include_icgc_data whether or not you want ICGC and other external data included
#' @param human_friendly_name a slightly more verbose name for your project
#' @param project_name unique ID for your project
#' @param description A verbose description of your data set
#' @param out_dir The full path to the base directory where the files are being created
#' @param overwrite Flag to specify that files should be overwritten if they exist
#'
#' @return a vector of sample_id for the patients that have been included
#' @export
#'
#' @examples
setup_study = function(short_name = "GAMBL",
                       include_icgc_data = FALSE,
                       human_friendly_name = "GAMBL data",
                       project_name = "gambl_minus_icgc",
                       description = "GAMBL data without ICGC",
                       cancer_type = "mixed",
                       out_dir,
                       overwrite = FALSE,
                       exclude_cohorts = c("FFPE_Benchmarking")){
  meta_study = paste0(out_dir, "meta_study.txt")
  meta_cancer_type = paste0(out_dir, "meta_cancer_type.txt")
  data_cancer_type = paste0(out_dir, "data_cancer_type.txt")
  meta_clinical_samples = paste0(out_dir, "meta_clinical_samples.txt")
  meta_clinical_patients = paste0(out_dir, "meta_clinical_patient.txt")

  caselist_cna = paste0(out_dir, "case_lists/cases_cna.txt")

  detailed_sv_out = paste0(out_dir, "annotated_svs_with_detail.tsv")

  meta_cna = paste0(out_dir, "meta_CNA.txt")
  gistic_cn_file = paste0(out_dir, "gistic_data_CNA.txt")
  meta_mutations = paste0(out_dir, "meta_mutations_extended.txt")
  data_mutations = "data_mutations_extended.maf"
  data_mutations_full = paste0(out_dir, "data_mutations_extended.maf")

  #determine what table to query and what restrictions to use for the MAF data
  if(include_icgc_data){
    maf_table = gambl_icgc_maf
  }else{
    maf_table = gambl_maf
  }

  #set up the new directory if necessary
  if (!file.exists(out_dir)){
    dir.create(out_dir)
    dir.create(paste0(out_dir, "case_lists"))
  }

  meta_study_content = paste0("type_of_cancer: ", cancer_type, "\n",
                              "cancer_study_identifier: ", project_name, "\n",
                              "name: ", human_friendly_name, "\n",
                              "short_name: ", short_name, "\n",
                              "description: ", description, "\n",
                              "add_global_case_list: true\n")
  cat(meta_study_content, file = meta_study)

  meta_cancer_type_content = paste0("genetic_alteration_type: CANCER_TYPE\n",
                                    "datatype: CANCER_TYPE\n",
                                    "data_filename: data_cancer_type.txt\n")
  cat(meta_cancer_type_content, file = meta_cancer_type)
  cat("brca-es0	Breast Invasive Carcinoma	breast,breast invasive	HotPink	Breast\n", file = data_cancer_type)

  meta_clinical_samples_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                         "genetic_alteration_type: CLINICAL\n",
                                         "datatype: SAMPLE_ATTRIBUTES\n",
                                         "data_filename: data_clinical_samples.txt")
  cat(meta_clinical_samples_content, file = meta_clinical_samples)

  meta_clinical_patients_content = paste0("cancer_study_identifier: ", project_name, "\n",
                                          "genetic_alteration_type: CLINICAL\n",
                                          "datatype: PATIENT_ATTRIBUTES\n",
                                          "data_filename: data_clinical_patient.txt")
  cat(meta_clinical_patients_content, file = meta_clinical_patients)

  meta_mutations_content = paste0("cancer_study_identifier: ", project_name,"\n",
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
    coding_ssms = get_coding_ssm(table_name = gambl_maf, exclude_cohort = exclude_cohorts)
    write_tsv(coding_ssms, data_mutations_full, na = "")
  }else{
    #read in the MAF instead
    coding_ssms = data.table::fread(
      file = data_mutations_full,
      sep = "\t",
      stringsAsFactors = FALSE,
      verbose = FALSE,
      data.table = TRUE,
      showProgress = TRUE,
      header = TRUE,
      fill = TRUE,
      skip = "Hugo_Symbol",
      quote = ""
    )
  }
  ids = coding_ssms %>% 
    pull(Tumor_Sample_Barcode) %>% 
    unique()
  return(ids)
}
