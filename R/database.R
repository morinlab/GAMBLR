
#' Get GAMBL metadata
#'
#' @param db The GAMBL database name
#' @param seq_type_filter Filtering criteria (default: all genomes)
#' @param tissue_status_filter Filtering criteria (default: only tumour genomes)
#' @param case_set optional short name for a pre-defined set of cases avoiding any
#' embargoed cases (current options: 'BLGSP-study', 'FL-DLBCL-study')
#'
#' @return A data frame with metadata for each biopsy in GAMBL
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
#' # basic usage
#' my_metadata = get_gambl_metadata()
#' # use pre-defined custom sample sets
#' only_blgsp_metadata = get_gambl_metadata(case_set="BLGSP-study")
#' # override default filters and request metadata for samples other than tumour genomes, e.g. also get the normals
#' only_normal_metadata = get_gambl_metadata(tissue_status_filter = c('tumour','normal'))
get_gambl_metadata = function(db="gambl_test",seq_type_filter = "genome",
                              tissue_status_filter=c("tumour"), case_set){
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  sample_meta = tbl(con,"sample_metadata") %>% filter(seq_type == seq_type_filter & tissue_status %in% tissue_status_filter)
  #if we only care about genomes, we can drop/filter anything that isn't a tumour genome
  #The key for joining this table to the mutation information is to use sample_id. Think of this as equivalent to a library_id. It will differ depending on what assay was done to the sample.
  biopsy_meta = tbl(con,"biopsy_metadata") %>% select(-patient_id) %>% select(-pathology) %>% select(-time_point) %>% select(-EBV_status_inf) #drop duplicated columns
  all_meta = left_join(sample_meta,biopsy_meta,by="biopsy_id") %>% as.data.frame()
  all_meta[all_meta$pathology=="B-cell unclassified","pathology"] = "HGBL"  #TODO fix this in the metadata
  if(!missing(case_set)){
    if(case_set == "FL-DLBCL-study"){
      #get FL cases and DLBCL cases not in special/embargoed cohorts
      all_meta = all_meta %>% filter(pathology %in% c("FL","DLBCL")) %>% filter(!cohort %in% c("DLBCL_ctDNA","DLBCL_BLGSP","LLMPP_P01","DLBCL_LSARP_Trios"))
    }
    if(case_set == "BLGSP-study"){
      #get BL cases minus duplicates (i.e. drop benchmarking cases)
      all_meta = all_meta %>% filter(cohort %in% c("BL_Adult","BL_cell_lines","BL_ICGC","BLGSP_Bcell_UNC","BL_Pediatric"))
    }else if(case_set == "GAMBL-all"){
      #get all GAMBL but remove FFPE benchmarking cases and ctDNA
      all_meta = all_meta %>% filter(!cohort %in% c("FFPE_Benchmarking","DLBCL_ctDNA"))
    }
  }
  #add some derivative columns that simplify and consolidate some of the others (DLBCL-specific)
  all_meta = all_meta %>% mutate(lymphgen = case_when(
    pathology != "DLBCL" ~ pathology,
    str_detect(lymphgen_cnv_noA53,"/") ~ "COMPOSITE",
    TRUE ~ lymphgen_cnv_noA53
  ))


  all_meta = all_meta %>% mutate(consensus_coo_dhitsig = case_when(
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
  all_meta = all_meta %>% mutate(pathology_rank = case_when(
    pathology == "B-ALL" ~ 0,
    pathology == "CLL" ~ 1,
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
  all_meta = all_meta %>% mutate(lymphgen_rank = case_when(
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
  DBI::dbDisconnect(con)
  return(all_meta)
}

#' Get the patient-centric clinical metadata
#'
#' @param db optional. Specify a different database (changing this is not recommended)
#' @param time_unit Return follow-up times in one of three time units: year, month or day
#' @param censor_cbioportal Optionally request the censoring to be encoded in the specific style required by cBioPortal
#'
#' @return Data frame with one row for each patient_id
#' @export
#' @import tidyverse RMariaDB DBI
#'
#' @examples
#' outcome_df = get_gambl_outcomes()
get_gambl_outcomes = function(db="gambl_test",patient_ids,time_unit="year",censor_cbioportal=FALSE,complete_missing=FALSE){
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  all_outcome = tbl(con,"outcome_metadata") %>% as.data.frame()
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
#' @param db The GAMBL database name
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
#' @import DBI RMariaDB tidyverse
#'
#' @examples
#' #lazily get every SV in the table with default quality filters
#' all_sv = get_manta_sv()
#' #get all SVs for a single sample
#' some_sv = get_manta_sv(sample_id="94-15772_tumorA")
#' #get the SVs in a region around MYC
#' myc_locus_sv = get_manta_sv(region="8:128723128-128774067")
get_manta_sv = function(db="gambl_test",table_name="bedpe_manta_hg19",min_vaf=0.1,min_score=40,pass=TRUE,pairing_status,sample_id,chromosome,qstart,qend,region,with_chr_prefix=FALSE){

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
    all_sv = tbl(con,table_name) %>%
      filter((CHROM_A == chromosome & START_A >= qstart & START_A <= qend) | (CHROM_B == chromosome & START_B >= qstart & START_B <= qend)) %>%
      filter(VAF_tumour >= min_vaf & SOMATIC_SCORE >= min_score)
  }else{
    all_sv = tbl(con,table_name) %>% filter(VAF_tumour >= min_vaf & SOMATIC_SCORE >= min_score)
  }
  if(pass){
    all_sv = all_sv %>% filter(FILTER == "PASS")
  }
  if(!missing(pairing_status)){
    all_sv = all_sv %>% filter(pair_status == pairing_status)
  }
  if(!missing(sample_id)){
    all_sv = all_sv %>% filter(tumour_sample_id == sample_id)
  }
  all_sv = as.data.frame(all_sv)
  if(with_chr_prefix){
    #add chr prefix only if it's missing

    all_sv = all_sv %>% mutate(CHROM_A = case_when(
      str_detect(CHROM_A,"chr") ~ CHROM_A,
      TRUE ~ paste0("chr",CHROM_A)
    ))
    all_sv = all_sv %>% mutate(CHROM_B = case_when(
      str_detect(CHROM_B,"chr") ~ CHROM_B,
      TRUE ~ paste0("chr",CHROM_B)
    ))

  }
  DBI::dbDisconnect(con)
  return(all_sv)
}

#' Retrieve the nearest SV to a specified SV in a given patient
#'
#' @param db The GAMBL database name
#' @param table_name The table we are querying
#' @param min_vaf The minimum tumour VAF for a SV to be returned
#' @param min_score The lowest Manta somatic score for a SV to be returned
#' @param pair_status Use to restrict results (if desired) to matched or unmatched results (default is to return all)
#' @param sv_data A single row of a data frame returned by a function such as get_manta_sv
#' @param with_chr_prefix Prepend all chromosome names with chr (required by some downstream analyses)
#'
#' @return A data frame in a bedpe-like format with additional columns that allow filtering of high-confidence SVs
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
fetch_next_sv = function(db="gambl_test",table_name="bedpe_manta_hg19",min_vaf=0.1,min_score=40,pass=TRUE,sv_data,with_chr_prefix=FALSE,in_from){
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  #in from B, out from A
  if(in_from == "A"){
    #search for the nearest SV in the positive direction if the strand of the breakpoint is +, otherwise search in the negative direction.
    #Logic is reversed if from == A
    chromosome = sv_data$CHROM_B #searching from the "right" end of the SV

    sample = sv_data$tumour_sample_id
    all_sv = tbl(con,table_name) %>% filter(tumour_sample_id == sample) %>% filter(VAF_tumour >= min_vaf & SOMATIC_SCORE >= min_score) %>% as.data.frame()
    DBI::dbDisconnect(con)

    if(sv_data$STRAND_B=="+"){
      #get the next nearest SV with coordinate greater than this one
      qstart = sv_data$START_B
      all_sv_A = all_sv %>% filter(CHROM_A == chromosome & START_A > qstart & STRAND_A == "+") %>% arrange(START_A) %>% head(1)
      all_sv_B = all_sv %>% filter(CHROM_B == chromosome & START_B > qstart & STRAND_B == "+") %>% arrange(START_B) %>% head(1)
      print(all_sv_A)
      print(all_sv_B)
      #the only strand compatible with this is the same strand for A

    }else{
      qstart = sv_data$START_B
      all_sv_A = all_sv %>% filter(CHROM_A == chromosome & START_A < qstart & STRAND_A == "-") %>% arrange(START_A) %>% tail(1)
      all_sv_B = all_sv %>% filter(CHROM_B == chromosome & START_B < qstart & STRAND_B == "-") %>% arrange(START_B) %>% tail(1)
      print(all_sv_A)
      print(all_sv_B)
      closest = which.min(c(all_sv_A$START_A,all_sv_B$START_B))
    }
  }else{ #in_from == B
    chromosome = sv_data$CHROM_A #searching from the "left" end of the SV
    if(sv_data$STRAND_A=="+"){
      #get the next nearest SV with coordinate greater than this one
      qstart = sv_data$START_A
      all_sv_A = all_sv %>% filter(CHROM_A == chromosome & START_A > qstart & STRAND_A == "+") %>% arrange(START_A) %>% head(1)
      all_sv_B = all_sv %>% filter(CHROM_B == chromosome & START_B > qstart & STRAND_B == "+") %>% arrange(START_B) %>% head(1)
      print(all_sv_A)
      print(all_sv_B)
      closest = which.min(c(all_sv_A$START_A,all_sv_B$START_B))
    }else{
      qstart = sv_data$START_A
      all_sv_A = all_sv %>% filter(CHROM_A == chromosome & START_A < qstart & STRAND_A == "-") %>% arrange(START_A) %>% tail(1)
      all_sv_B = all_sv %>% filter(CHROM_B == chromosome & START_B < qstart & STRAND_B == "-") %>% arrange(START_B) %>% tail(1)
      print(all_sv_A)
      print(all_sv_B)
      closest = which.min(c(all_sv_A$START_A,all_sv_B$START_B))
    }
  }
  if(closest == 1){
    print("closest was 1")
    return(list(sv=all_sv_A,other_end="B"))
  }else{
    print("closest was 2")
    return(list(sv=all_sv_B,other_end="A"))
  }

}

#' Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range
#'
#' @param db The GAMBL database name
#' @param table_name The table we are querying
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
#' prefixed_segments = get_cn_segments(get_cn_segments(chromosome="12",qstart=122456912,qend=122464036,with_chr_prefix = TRUE))
get_cn_segments = function(db="gambl_test",table_name="seg_battenberg_hg19",chromosome,qstart,qend,region,with_chr_prefix=FALSE){
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
  all_segs = tbl(con,table_name) %>%
    filter((chrom == chromosome & start <= qstart & end >= qend) |
             (chrom == chromosome & start >= qstart & end <= qend)) %>%
    as.data.frame()
  if(! with_chr_prefix){
    all_segs = all_segs %>% mutate(chrom = gsub("chr","",chrom))
  }
  DBI::dbDisconnect(con)
  return(all_segs)
}



#' Get all somatic mutations for a given gene or list of genes and optionally restrict to coding variants
#'
#' @param db The GAMBL database name
#' @param table_name The table you are querying
#' @param gene_symbol Character vector of gene symbols
#'
#' @return MAF-format data frame of mutations in query gene
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
#' #basic usage
#' get_ssm_by_gene(gene_symbol=c("EZH2"),coding_only=TRUE)
get_ssm_by_gene = function(db="gambl_test",table_name="maf_slms3_hg19_icgc",gene_symbol,coding_only=FALSE){
  coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  muts_gene = tbl(con,table_name) %>%
    filter(Hugo_Symbol %in% gene_symbol)
  if(coding_only){
    muts_gene = muts_gene %>% filter(Variant_Classification %in% coding_class)
  }
  muts_gene = as.data.frame(muts_gene)
  DBI::dbDisconnect(con)
  return(muts_gene)
}

#' Retrieve all SSMs from the GAMBL database within a single genomic coordinate range
#'
#' @param db The GAMBL database name
#' @param table_name The table we are querying
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix)
#' @param qstart Query start coordinate of the range you are restricting to
#' @param qend Query end coordinate of the range you are restricting to
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately
#' @param basic_columns Set to TRUE to override the default behaviour of returning only the first 45 columns of MAF data
#'
#' @return A data frame containing all the MAF data columns (one row per mutation)
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
#' #basic usage
#' my_mutations=get_ssm_by_region(region="chr8:128,723,128-128,774,067")
#' #specifying chromosome, start and end individually
#' my_mutations=get_ssm_by_region(chromosome="8",qstart=128723128,qend=128774067)
get_ssm_by_region = function(db="gambl_test",table_name="maf_slms3_hg19_icgc",chromosome,qstart,qend,region="",basic_columns=TRUE,streamlined){
  if(!region==""){
    region = gsub(",","",region)
    #format is chr6:37060224-37151701
    split_chunks = unlist(strsplit(region,":"))
    chromosome = split_chunks[1]

    startend = unlist(strsplit(split_chunks[2],"-"))
    qstart=startend[1]
    qend=startend[2]
  }
  chromosome = gsub("chr","",chromosome)

  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  muts_region = tbl(con,table_name) %>%
    filter(Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)

  muts_region = as.data.frame(muts_region)
  if(!missing(streamlined)){
    muts_region = muts_region %>% select(Start_Position,Tumor_Sample_Barcode)
  }else if(basic_columns){
    muts_region = muts_region[,c(1:45)]
  }
  DBI::dbDisconnect(con)
  return(muts_region)
}

#' Retrieve all coding SSMs from the GAMBL database in MAF-like format
#'
#' @param db The GAMBL database name
#' @param table_name The table we are querying
#' @param limit_cohort Supply this to restrict mutations to one or more cohorts in a list
#' @param exclude_cohort  Supply this to exclude mutations from one or more cohorts in a list
#' @param limit_pathology Supply this to restrict mutations to one pathology
#' @param basic_columns Set to TRUE to override the default behaviour of returning only the first 45 columns of MAF data
#'
#' @return A data frame containing all the MAF data columns (one row per mutation)
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
#' #basic usage
#' maf_data = get_coding_ssm(limit_cohort=c("BL_Adult","BL_Pediatric","BL_ICGC"))
get_coding_ssm = function(db="gambl_test",table_name="maf_slms3_hg19_icgc",limit_cohort,exclude_cohort,limit_pathology,basic_columns=TRUE){
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")
  sample_meta = tbl(con,"sample_metadata") %>% filter(seq_type == "genome" & tissue_status == "tumour")
  biopsy_meta = tbl(con,"biopsy_metadata") %>% select(-patient_id) %>%
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
  muts = tbl(con,table_name) %>%
    filter(Variant_Classification %in% coding_class & Tumor_Sample_Barcode %in% sample_ids) %>%
   as.data.frame()

  if(basic_columns){
    muts = muts[,c(1:45)]
  }
  DBI::dbDisconnect(con)
  return(muts)
}


#TODO migrate viz/plotting functions that don't directly rely on the database to a separate file
#' Make a rainbow plot of all mutations in a region, ordered and coloured by metadata
#'
#' @param db The GAMBL database name
#' @param table_name The table we are querying
#' @param mutations_maf A data frame containing mutations (MAF format) within a region of interest (i.e. use the get_ssm_by_region)
#' @param metadata should be a data frame with sample_id as a column that should match Tumor_Sample_Barcode in the database
#' @param classification_column The name of the metadata column to use for ordering and colouring samples
#' @param bed Optional data frame specifying the regions to annotate (required columns: start, end, name)
#'
#' @return ggplot2 object
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
#' # basic usage
#' ashm_rainbow_plot(mutations_maf=my_mutations,metadata=my_metadata)
#' # advanced usages
#' mybed = data.frame(start=c(128806578,128805652,128748315), end=c(128806992,128809822,128748880), name=c("TSS","enhancer","MYC-e1"))
#' ashm_rainbow_plot(mutations_maf=my_mutations,metadata=my_metadata,bed=mybed)
ashm_rainbow_plot = function(db="gambl_test",table_name="maf_slms3_hg19_icgc",mutations_maf,metadata,classification_column="pathology",bed){
  mutation_positions = mutations_maf %>%
    select(Tumor_Sample_Barcode,Start_Position) %>% as.data.frame()
  muts_anno = left_join(mutation_positions,metadata,by=c("Tumor_Sample_Barcode" = "sample_id"))
  ordering = metadata$sample_id[order(metadata[,classification_column])]
  muts_anno$sample_id = factor(muts_anno$Tumor_Sample_Barcode,levels=ordering)
  muts_anno$classification = muts_anno[,classification_column]

  p = ggplot(muts_anno) + geom_point(aes(x=Start_Position,y=sample_id,colour=classification),alpha=0.4) + theme(axis.text.y=element_blank())
  if(missing(bed)){
    p
  }else{
    height = length(unique(muts_anno$Tumor_Sample_Barcode))
    p = p + geom_rect(data=bed, aes(xmin = start, xmax = end, ymin = 0, ymax = height+30),alpha=0.1) + geom_text(data=bed,aes(x = start, y= height,label=name),size = 3,angle=90)
  }
  return(p)
}

#' Title
#'
#' @param regions_bed Bed file with chromosome coordinates, should contain columns chr, start, end, name (with these exact names)
#' @param regions_to_display Optional vector of names from default regions_bed to use
#' @param metadata A metadata file already subsetted and arranged on the order you want the samples vertically displayed
#' @param classification_column optional. Override default column for assigning the labels used for colouring in the figure.
#' @param db optional. Specify a different database (changing this is not recommended)
#' @param table_name optional. Specify a different table to query (changing this is not recommended)
#'
#' @return nothing
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
ashm_multi_rainbow_plot = function(regions_bed,regions_to_display,exclude_classifications,metadata,custom_colours,classification_column="lymphgen",db="gambl_test",table_name="maf_slms3_hg19_icgc"){
  #get the mutations for each region and combine
  #regions_bed should contain chr, start, end, name (with these exact names)
  if(missing(metadata)){
    metadata = get_gambl_metadata()
    meta_arranged = arrange(metadata,pathology_rank,lymphgen)
  }else{
    meta_arranged = metadata #assume the user already arranged it the way they wanted
  }
  if(!missing(exclude_classifications)){
    meta_arranged = filter(meta_arranged,!get(classification_column) %in% exclude_classifications)
  }
  if(missing(regions_bed)){
    regions_bed= grch37_ashm_regions
    regions_bed = mutate(regions_bed,regions=paste0(chr_name,":",hg19_start,"-",hg19_end))
    regions_bed = mutate(regions_bed,name=paste0(gene,"-",region))
  }else{
    regions_bed = mutate(regions_bed,regions=paste0(chr,":",start,"-",end))
  }

  names=pull(regions_bed,name)
  names = c(names,"CCND1","FOXP1-TSS1","FOXP1-TSS2","FOXP1-TSS3","FOXP1-TSS4","FOXP1-TSS5","BCL6","IGH","IGL","IGK","PVT1","BCL2") #add some additional regions of interest
  regions = pull(regions_bed,regions)
  regions = c(regions,"chr11:69451233-69460334","chr3:71623481-71641671","chr3:71532613-71559445","chr3:71343345-71363145","chr3:71167050-71193679","chr3:71105715-71118362",
              "chr3:187406804-188522799","chr14:106144562-106344765","chr22:23217074-23250428","chr2:89073691-89320640","chr8:128774985-128876311","chr18:60982124-60990180")

  region_mafs = lapply(regions,function(x){get_ssm_by_region(region=x,streamlined = TRUE)})
  tibbled_data = tibble(region_mafs, region_name = names)
  unnested_df = tibbled_data %>% unnest_longer(region_mafs)
  unlisted_df = mutate(unnested_df,start=region_mafs$Start_Position,sample_id=region_mafs$Tumor_Sample_Barcode) %>%
      select(start,sample_id,region_name)



  meta_arranged$classification = factor(meta_arranged[,classification_column],levels=unique(meta_arranged[,classification_column]))
  muts_anno = left_join(unlisted_df,meta_arranged)
  muts_first =  select(muts_anno,start,region_name) %>% group_by(region_name) %>% arrange(start) %>% filter(row_number()==1)
  eg = expand_grid(start=pull(muts_first,start),sample_id=pull(meta_arranged,sample_id))
  eg = left_join(eg,muts_first)

  #concatenate expanded frame of points with original mutation data
  real_and_fake = bind_rows(unlisted_df,eg)
  muts_anno = left_join(real_and_fake,meta_arranged) #%>% filter(!is.na(get(classification_column)))

  muts_anno$sample_id= factor(muts_anno$sample_id,levels=meta_arranged$sample_id)

  if(!missing(regions_to_display)){
    muts_anno = filter(muts_anno,region_name %in% regions_to_display)
  }
  #make the plot
  if(missing(custom_colours)){
    p = muts_anno %>%
      ggplot() + geom_point(aes(x=start,y=sample_id,colour=classification),alpha=0.4,size=0.6) +
                    theme(axis.text.y=element_blank()) +
    facet_wrap(~region_name,scales="free_x") + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))

    print(p)
  }else{
  #testing manual colouring
  p = muts_anno %>%
    ggplot() + geom_point(aes(x=start,y=sample_id,colour=classification),alpha=0.4,size=0.6) +
    theme(axis.text.y=element_blank()) + scale_colour_manual(values=custom_colours) +
    facet_wrap(~region_name,scales="free_x") + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
  print(p)
  }
}

#' Create a genome-wide copy number plot for one sample and (optionally) display mutation VAF
#'
#' @param this_sample
#' @param just_segments
#' @param genes_to_label optional. Provide a list of genes to label (if mutated). Can only be used with coding_only (see below)
#' @param coding_only optional. Set to TRUE to restrict to plotting only coding mutations
#'
#' @return nothing
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
copy_number_vaf_plot = function(this_sample,just_segments=FALSE,coding_only=FALSE,genes_to_label){
  chrom_order=factor(c(1:22,"X"))
  cn_colours = get_gambl_colours(classification = "copy_number")
  maf_and_seg = assign_cn_to_ssm(this_sample=this_sample,coding_only=coding_only)
  vaf_cn_maf = maf_and_seg[["maf"]]
  vaf_cn_maf = mutate(vaf_cn_maf,CN=as.character(CN))
  if(just_segments){
    #I realized this is ugly
    cn_seg = maf_and_seg[["seg"]]
    cn_seg = mutate(cn_seg,CN_segment = as.numeric(CN),CN = as.character(CN))
    ggplot(cn_seg) +
      geom_segment(data=cn_seg,aes(x=Start_Position,xend=End_Position,y=CN_segment,yend=CN_segment)) +
      facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") +
      scale_colour_manual(values = cn_colours) +
      theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
  }else{
    if(coding_only){
      if(missing(genes_to_label)){
        p = mutate(vaf_cn_maf,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
          geom_point(aes(x=Start_Position,y=vaf,colour=CN),alpha=0.6,size=2) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") +
          theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
        p
      }else{
        #label any mutations that intersect with our gene list
        plot_genes = vaf_cn_maf %>% filter(Hugo_Symbol %in% my_genes)

        p = mutate(vaf_cn_maf,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
          geom_point(aes(x=Start_Position,y=vaf,colour=CN),size=2) +
          geom_text(data=plot_genes,aes(x=Start_Position,y=0.8,label=Hugo_Symbol),size=3,angle=90) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") + ylim(c(0,1)) +
          theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
        p
      }
    }else{

      p = mutate(vaf_cn_maf,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
         geom_point(aes(x=Start_Position,y=vaf,colour=CN),alpha=0.6,size=0.2) +
         scale_colour_manual(values = cn_colours) +
         facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") +
         theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
      p
    }
  }
}


