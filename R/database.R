
#' Get GAMBL metadata
#'
#' @param seq_type_filter Filtering criteria (default: all genomes)
#' @param tissue_status_filter Filtering criteria (default: only tumour genomes)
#' @param case_set optional short name for a pre-defined set of cases avoiding any
#' embargoed cases (current options: 'BLGSP-study', 'FL-DLBCL-study', 'DLBCL-unembargoed)
#'
#' @return A data frame with metadata for each biopsy in GAMBL
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' # basic usage
#' my_metadata = get_gambl_metadata()
#' # use pre-defined custom sample sets
#' only_blgsp_metadata = get_gambl_metadata(case_set="BLGSP-study")
#' # override default filters and request metadata for samples other than tumour genomes, e.g. also get the normals
#' only_normal_metadata = get_gambl_metadata(tissue_status_filter = c('tumour','normal'))
get_gambl_metadata = function(seq_type_filter = "genome",
                              tissue_status_filter=c("tumour"), case_set, remove_benchmarking = TRUE, with_outcomes=FALSE){
  db=config::get("database_name")
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  sample_meta = dplyr::tbl(con,"sample_metadata")
  sample_meta_normal_genomes =  sample_meta %>% filter(seq_type == "genome" & tissue_status=="normal") %>%
    select(patient_id,sample_id) %>% as.data.frame() %>% rename("normal_sample_id"="sample_id")

  sample_meta = sample_meta %>% filter(seq_type == seq_type_filter & tissue_status %in% tissue_status_filter)

  #if we only care about genomes, we can drop/filter anything that isn't a tumour genome
  #The key for joining this table to the mutation information is to use sample_id. Think of this as equivalent to a library_id. It will differ depending on what assay was done to the sample.
  biopsy_meta = dplyr::tbl(con,"biopsy_metadata") %>% select(-patient_id) %>% select(-pathology) %>% select(-time_point) %>% select(-EBV_status_inf) #drop duplicated columns
  all_meta = dplyr::left_join(sample_meta,biopsy_meta,by="biopsy_id") %>% as.data.frame()
  if(seq_type_filter == "genome" & length(tissue_status_filter) == 1 & tissue_status_filter[1] == "tumour"){
    #join back the matched normal genome
    all_meta = left_join(all_meta,sample_meta_normal_genomes,by="patient_id")
    all_meta = all_meta %>% mutate(pairing_status=case_when(is.na(normal_sample_id)~"unmatched",TRUE~"matched"))
  }


  #all_meta[all_meta$pathology=="B-cell unclassified","pathology"] = "HGBL"  #TODO fix this in the metadata
  if(remove_benchmarking){
    all_meta = all_meta %>% filter(cohort != "FFPE_Benchmarking")
  }
  if(!missing(case_set)){
    if(case_set == "FL-DLBCL-study"){
      #get FL cases and DLBCL cases not in special/embargoed cohorts
      all_meta = all_meta %>% dplyr::filter(pathology %in% c("FL","DLBCL")) %>% filter(!cohort %in% c("DLBCL_ctDNA","DLBCL_BLGSP","LLMPP_P01","DLBCL_LSARP_Trios"))
    }
    if(case_set == "DLBCL-unembargoed"){
      #get DLBCL cases not in special/embargoed cohorts
      all_meta = all_meta %>% dplyr::filter(pathology %in% c("DLBCL")) %>% filter(!cohort %in% c("DLBCL_ctDNA","DLBCL_BLGSP","LLMPP_P01","DLBCL_LSARP_Trios","DLBCL_HTMCP"))
    }
    if(case_set == "BLGSP-study"){
      #get BL cases minus duplicates (i.e. drop benchmarking cases)
      all_meta = all_meta %>% dplyr::filter(cohort %in% c("BL_Adult","BL_cell_lines","BL_ICGC","BLGSP_Bcell_UNC","BL_Pediatric") |(cohort=="LLMPP_P01" & pathology == "BL"))
    }else if(case_set == "GAMBL-all"){
      #get all GAMBL but remove FFPE benchmarking cases and ctDNA
      all_meta = all_meta %>% dplyr::filter(!cohort %in% c("FFPE_Benchmarking","DLBCL_ctDNA"))
    }
  }



  #add some derivative columns that simplify and consolidate some of the others (DLBCL-specific)
  all_meta = all_meta %>% dplyr::mutate(lymphgen = case_when(
    pathology != "DLBCL" & pathology != "FL" ~ pathology,
    str_detect(lymphgen_cnv_noA53,"/") ~ "COMPOSITE",
    TRUE ~ lymphgen_cnv_noA53
  ))

  all_meta = mutate(all_meta,Tumor_Sample_Barcode=sample_id) #duplicate for convenience
  all_meta = all_meta %>% dplyr::mutate(consensus_coo_dhitsig = case_when(
    pathology != "DLBCL" ~ pathology,
    COO_consensus == "ABC" ~ COO_consensus,
    DLBCL90_dhitsig_call == "POS" ~ "DHITsigPos",
    DLBCL90_dhitsig_call == "NEG" ~ "DHITsigNeg",
    DHITsig_PRPS_class == "DHITsigPos" ~ "DHITsigPos",
    DHITsig_PRPS_class == "DHITsigNeg" ~ "DHITsigNeg",
    DHITsig_PRPS_class == "UNCLASS" ~ "DHITsigPos",
    TRUE ~ "NA"
  ))
  #assign a rank to each pathology for consistent and sensible ordering
  all_meta = all_meta %>% dplyr::mutate(pathology_rank = case_when(
    pathology == "B-ALL" ~ 0,
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
    TRUE ~ 35
  ))
  all_meta = all_meta %>% dplyr::mutate(lymphgen_rank = case_when(
    pathology != "DLBCL" ~ pathology_rank,
    lymphgen == "Other" ~ 16,
    lymphgen == "COMPOSITE" ~ 17,
    lymphgen == "N1" ~ 18,
    lymphgen == "EZB" ~ 19,
    lymphgen == "ST2" ~ 20,
    lymphgen == "BN2" ~ 21,
    lymphgen == "MCD" ~ 22,
    TRUE ~ 50
  ))
  if(with_outcomes){
    outcome_table = get_gambl_outcomes() %>% select(-sex)
    all_meta = left_join(all_meta,outcome_table,by="patient_id") %>%
      mutate(age_group = case_when(cohort=="BL_Adult"~"Adult_BL",cohort=="BL_Pediatric" | cohort == "BL_ICGC" ~ "BL_Pediatric", TRUE ~ "Other"))

  }
  DBI::dbDisconnect(con)
  return(all_meta)
}

#' Get the patient-centric clinical metadata
#'
#' @param time_unit Return follow-up times in one of three time units: year, month or day
#' @param censor_cbioportal Optionally request the censoring to be encoded in the specific style required by cBioPortal
#'
#' @return Data frame with one row for each patient_id
#' @export
#' @import tidyverse RMariaDB DBI dbplyr
#'
#' @examples
#' outcome_df = get_gambl_outcomes()
get_gambl_outcomes = function(patient_ids,time_unit="year",censor_cbioportal=FALSE,complete_missing=FALSE){
  db=config::get("database_name")
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  all_outcome = dplyr::tbl(con,"outcome_metadata") %>% as.data.frame()
  if(!missing(patient_ids)){
    all_outcome = all_outcome %>% filter(patient_id %in% patient_ids)
    if(complete_missing){
      #add NA values and censored outcomes for all missing patient_ids
      all_outcome = all_outcome %>% complete(patient_id= patient_ids,fill=list(OS_YEARS=0,PFS_years=0,TTP_YEARS=0,DSS_YEARS=0,CODE_OS=0,CODE_PFS=0,CODE_DSS=0,CODE_TTP=0))
    }
  }
  if(time_unit == "month"){
    all_outcome = all_outcome %>% mutate(OS_MONTHS=OS_YEARS * 12)
    all_outcome = all_outcome %>% mutate(PFS_MONTHS=PFS_YEARS * 12)
    all_outcome = all_outcome %>% mutate(TTP_MONTHS=TTP_YEARS * 12)
    all_outcome = all_outcome %>% mutate(DSS_MONTHS=DSS_YEARS * 12)
    all_outcome = all_outcome %>% select(-c("OS_YEARS","PFS_YEARS","TTP_YEARS","DSS_YEARS"))
  }else if(time_unit == "day"){
    all_outcome = all_outcome %>% mutate(OS_DAYS=OS_YEARS * 365)
    all_outcome = all_outcome %>% mutate(PFS_DAYS=PFS_YEARS * 365)
    all_outcome = all_outcome %>% mutate(TTP_DAYS=TTP_YEARS * 365)
    all_outcome = all_outcome %>% mutate(DSS_DAYS=DSS_YEARS * 365)
    all_outcome = all_outcome %>% select(-c("OS_YEARS","PFS_YEARS","TTP_YEARS","DSS_YEARS"))
  }
  #if necessary, convert the censoring into the cBioPortal format for OS and PFS
  if(censor_cbioportal){
    all_outcome$OS_STATUS = as.character(all_outcome$CODE_OS)
    all_outcome = all_outcome %>% mutate(OS_STATUS = case_when(OS_STATUS=="0" ~ "0:LIVING",OS_STATUS=="1"~"1:DECEASED"))
    all_outcome$DFS_STATUS = as.character(all_outcome$CODE_PFS)
    all_outcome = all_outcome %>% mutate(DFS_STATUS = case_when(DFS_STATUS=="0" ~ "0:DiseaseFree",DFS_STATUS=="1"~"1:Recurred/Progressed"))
    all_outcome = all_outcome %>% mutate(all_outcome,DFS_MONTHS=PFS_MONTHS)
  }
  DBI::dbDisconnect(con)
  return(all_outcome)
}

#' Retrieve Manta SVs from the database and filter
#'
#' @param table_name The table we are querying
#' @param min_vaf The minimum tumour VAF for a SV to be returned
#' @param min_score The lowest Manta somatic score for a SV to be returned
#' @param pair_status Use to restrict results (if desired) to matched or unmatched results (default is to return all)
#' @param chromosome The chromosome you are restricting to
#' @param qstart Query start coordinate of the range you are restricting to
#' @param qend Query end coordinate of the range you are restricting to
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately
#' @param with_chr_prefix Prepend all chromosome names with chr (required by some downstream analyses)
#'
#' @return A data frame in a bedpe-like format with additional columns that allow filtering of high-confidence SVs
#' @export
#' @import DBI RMariaDB tidyverse dbplyr
#'
#' @examples
#' #lazily get every SV in the table with default quality filters
#' all_sv = get_manta_sv()
#' #get all SVs for a single sample
#' some_sv = get_manta_sv(sample_id="94-15772_tumorA")
#' #get the SVs in a region around MYC
#' myc_locus_sv = get_manta_sv(region="8:128723128-128774067")
get_manta_sv = function(min_vaf=0.1,min_score=40,pass=TRUE,pairing_status,sample_id,chromosome,qstart,qend,region,with_chr_prefix=FALSE){
  db=config::get("database_name")
  table_name=config::get("results_tables")$sv
    if(!missing(region)){
      region = gsub(",","",region)
      #format is chr6:37060224-37151701
      split_chunks = unlist(strsplit(region,":"))
      chromosome = split_chunks[1]
      startend = unlist(strsplit(split_chunks[2],"-"))
      qstart=startend[1]
      qend=startend[2]
    }
  #this table stores chromosomes with un-prefixed names. Convert to prefixed chromosome if necessary

  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  if(!missing(region) || !missing(chromosome)){
    if(grepl("chr",chromosome)){
      chromosome = gsub("chr","",chromosome)
    }
    all_sv = dplyr::tbl(con,table_name) %>%
      dplyr::filter((CHROM_A == chromosome & START_A >= qstart & START_A <= qend) | (CHROM_B == chromosome & START_B >= qstart & START_B <= qend)) %>%
      dplyr::filter(VAF_tumour >= min_vaf & SOMATIC_SCORE >= min_score)
  }else{
    all_sv = dplyr::tbl(con,table_name) %>% dplyr::filter(VAF_tumour >= min_vaf & SOMATIC_SCORE >= min_score)
  }
  if(pass){
    all_sv = all_sv %>% dplyr::filter(FILTER == "PASS")
  }
  if(!missing(pairing_status)){
    all_sv = all_sv %>% dplyr::filter(pair_status == pairing_status)
  }
  if(!missing(sample_id)){
    all_sv = all_sv %>% dplyr::filter(tumour_sample_id == sample_id)
  }
  all_sv = as.data.frame(all_sv)
  if(with_chr_prefix){
    #add chr prefix only if it's missing

    all_sv = all_sv %>% dplyr::mutate(CHROM_A = case_when(
      str_detect(CHROM_A,"chr") ~ CHROM_A,
      TRUE ~ paste0("chr",CHROM_A)
    ))
    all_sv = all_sv %>% dplyr::mutate(CHROM_B = case_when(
      str_detect(CHROM_B,"chr") ~ CHROM_B,
      TRUE ~ paste0("chr",CHROM_B)
    ))

  }
  DBI::dbDisconnect(con)
  return(all_sv)
}


#' Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range
#'
#' @param chromosome The chromosome you are restricting to
#' @param start Start coordinate of the range you are restricting to
#' @param end End coordinate of the range you are restricting to
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately
#' @param with_chr_prefix Prepend all chromosome names with chr (required by some downstream analyses)
#'
#' @return A data frame containing all the MAF data columns (one row per mutation)
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
#' # basic usage
#' my_segments=get_cn_segments(region="chr8:128,723,128-128,774,067")
#' # specifying chromosome, start and end individually
#' my_segments=get_cn_segments(chromosome="8",qstart=128723128,qend=128774067)
#' # Asking for chromosome names to have a chr prefix (default is un-prefixed)
#' prefixed_segments = get_cn_segments(chromosome="12",qstart=122456912,qend=122464036,with_chr_prefix = TRUE)
get_cn_segments = function(chromosome,qstart,qend,region,with_chr_prefix=FALSE){
  db = config::get("database_name")
  table_name = config::get("results_tables")$copy_number
  if(!missing(region)){
    region = gsub(",","",region)
    #format is chr6:37060224-37151701
    split_chunks = unlist(strsplit(region,":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2],"-"))
    qstart=startend[1]
    qend=startend[2]
  }
  #chr prefix the query chromosome to match how it's stored in the table.
  # This isn't yet standardized in the db so it's just a workaround "for now".
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)

  if(grepl("chr",chromosome)){

  }else{
    chromosome = paste0("chr",chromosome)
  }
  #remove the prefix if this is false (or leave as is otherwise)

  #TODO improve this query to allow for partial overlaps
  all_segs = dplyr::tbl(con,table_name) %>%
    filter((chrom == chromosome & start <= qstart & end >= qend) |
             (chrom == chromosome & start >= qstart & end <= qend)) %>%
    as.data.frame()
  if(! with_chr_prefix){
    all_segs = all_segs %>% dplyr::mutate(chrom = gsub("chr","",chrom))
  }
  DBI::dbDisconnect(con)
  return(all_segs)
}

#' Title
#'
#' @param table_name
#' @param connection
#' @param file_path
#'
#' @return
#' @export
#'
#' @examples
append_to_table = function(table_name,data_df){
  db=config::get("database_name")
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  dbWriteTable(con,table_name,table_data,append=TRUE)
}


#' Get all somatic mutations for a given gene or list of genes and optionally restrict to coding variants
#'
#' @param gene_symbol Character vector of gene symbols
#'
#' @return MAF-format data frame of mutations in query gene
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' #basic usage
#' get_ssm_by_gene(gene_symbol=c("EZH2"),coding_only=TRUE)
get_ssm_by_gene = function(gene_symbol,coding_only=FALSE,rename_splice_region=TRUE){
  table_name = config::get("results_tables")$ssm
  db=config::get("database_name")
  coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  muts_gene = dplyr::tbl(con,table_name) %>%
    dplyr::filter(Hugo_Symbol %in% gene_symbol)
  if(coding_only){
    muts_gene = muts_gene %>% dplyr::filter(Variant_Classification %in% coding_class)
  }
  muts_gene = as.data.frame(muts_gene)
  if(rename_splice_region){
    muts_gene = muts_gene %>% mutate(Variant_Classification = case_when(
      Variant_Classification == "Splice_Region" ~ "Splice_Site",
      TRUE ~ Variant_Classification))
  }
  DBI::dbDisconnect(con)
  return(muts_gene)
}

#' Retrieve all SSMs from the GAMBL database within a single genomic coordinate range
#'
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix)
#' @param qstart Query start coordinate of the range you are restricting to
#' @param qend Query end coordinate of the range you are restricting to
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately
#' @param basic_columns Set to TRUE to override the default behaviour of returning only the first 45 columns of MAF data
#'
#' @return A data frame containing all the MAF data columns (one row per mutation)
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' #basic usage
#' my_mutations=get_ssm_by_region(region="chr8:128,723,128-128,774,067")
#' #specifying chromosome, start and end individually
#' my_mutations=get_ssm_by_region(chromosome="8",qstart=128723128,qend=128774067)
get_ssm_by_region = function(chromosome,qstart,qend,region="",basic_columns=TRUE,streamlined,maf_data){
  table_name = config::get("results_tables")$ssm
  db=config::get("database_name")
  if(!region==""){
    region = gsub(",","",region)
    #format is chr6:37060224-37151701
    split_chunks = unlist(strsplit(region,":"))
    chromosome = split_chunks[1]

    startend = unlist(strsplit(split_chunks[2],"-"))
    qstart=as.numeric(startend[1])
    qend=as.numeric(startend[2])
    #print(class(qstart))
    #print(class(qend))
  }
  chromosome = gsub("chr","",chromosome)
  if(missing(maf_data)){
    con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
    muts_region = dplyr::tbl(con,table_name) %>%
      dplyr::filter(Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)
    muts_region = as.data.frame(muts_region)
  }else{
    print(paste0("filtering based on Chromosome == ",chromosome," Start_Position >", qstart, "& Start_Position < ", qend))
    muts_region = dplyr::filter(maf_data,Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)
    test = muts_region[,c(1:7)]
    print(head(test))
    muts_region = dplyr::filter(maf_data,Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)
    test = muts_region[,c(1:7)]
    print(head(test))
  }

  if(!missing(streamlined)){
    muts_region = muts_region %>% dplyr::select(Start_Position,Tumor_Sample_Barcode)
  }else if(basic_columns){
    muts_region = muts_region[,c(1:45)]
  }
  if(missing(maf_data)){
    DBI::dbDisconnect(con)
  }
  return(muts_region)
}

#' Retrieve all coding SSMs from the GAMBL database in MAF-like format
#'
#' @param limit_cohort Supply this to restrict mutations to one or more cohorts in a list
#' @param exclude_cohort  Supply this to exclude mutations from one or more cohorts in a list
#' @param limit_pathology Supply this to restrict mutations to one pathology
#' @param basic_columns Set to TRUE to override the default behaviour of returning only the first 45 columns of MAF data
#'
#' @return A data frame containing all the MAF data columns (one row per mutation)
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' #basic usage
#' maf_data = get_coding_ssm(limit_cohort=c("BL_ICGC"))
get_coding_ssm = function(limit_cohort,exclude_cohort,limit_pathology,basic_columns=TRUE){
  table_name = config::get("results_tables")$ssm
  db=config::get("database_name")
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")
  sample_meta = dplyr::tbl(con,"sample_metadata") %>% filter(seq_type == "genome" & tissue_status == "tumour")
  biopsy_meta = dplyr::tbl(con,"biopsy_metadata") %>% select(-patient_id) %>%
    select(-pathology) %>% select(-time_point) %>% select(-EBV_status_inf) #drop duplicated columns
    all_meta = left_join(sample_meta,biopsy_meta,by="biopsy_id") %>%
    as.data.frame()

  #do all remaining filtering on the metadata then add the remaining sample_id to the query
  if(!missing(limit_cohort)){
    all_meta = all_meta %>% filter(cohort %in% limit_cohort)
  }
  if(!missing(exclude_cohort)){
    all_meta = all_meta %>% filter(!cohort %in% exclude_cohort)
  }
  if(!missing(limit_pathology)){
    all_meta = all_meta %>% filter(pathology %in% limit_pathology)
  }
  sample_ids = pull(all_meta,sample_id)
  muts = tbl(con,table_name) %>% filter(Variant_Classification %in% coding_class) %>% as.data.frame()
  muts = muts %>% filter(Tumor_Sample_Barcode %in% sample_ids)

  if(basic_columns){
    muts = muts[,c(1:45)]
  }
  DBI::dbDisconnect(con)
  return(muts)
}

#' Title
#'
#' @param hugo_symbols
#' @param tidy_expression_data
#' @param metadata
#' @param join_with
#'
#' @return
#' @export
#'
#' @examples
get_gene_expression = function(hugo_symbols,tidy_expression_data,metadata,join_with="mrna"){

  database_name = config::get("database_name")
  if(missing(metadata)){
    if(join_with=="mrna"){
      metadata = get_gambl_metadata(seq_type_filter = "mrna")
    }else{
      metadata = get_gambl_metadata()
    }
  }
  metadata = metadata %>% select(sample_id)
  if(missing(hugo_symbols)){
    print("ERROR: supply at least one gene symbol")
  }
  #load the tidy expression data from the database
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  tidy_expression_data = tbl(con,"expression_vst_hg38")

  gene_expression_df = tidy_expression_data %>%
    filter(Hugo_Symbol %in% hugo_symbols) %>% as.data.frame()


  if(join_with=="mrna"){
    #join to metadata
    gene_expression_df = select(gene_expression_df,-genome_sample_id,-biopsy_id)
    expression_wider = pivot_wider(gene_expression_df,names_from=Hugo_Symbol,values_from=expression)
    expression_wider = left_join(metadata,expression_wider,by=c("sample_id"="mrna_sample_id"))
  }else{
    expression_wider = select(gene_expression_df,-mrna_sample_id,-biopsy_id) %>%
      filter(genome_sample_id !="NA") %>%
      pivot_wider(names_from=Hugo_Symbol,values_from=expression)
    expression_wider = left_join(metadata,expression_wider,by=c("sample_id"="genome_sample_id"))
  }
  return(expression_wider)
}

