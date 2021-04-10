
require("dbplyr")
require("tidyverse")

#' Get GAMBL metadata
#'
#' @param db The GAMBL database name
#' @param seq_type_filter Filtering criteria (default: all genomes)
#' @param tissue_status_filter Filtering criteria (default: only tumour genomes)
#'
#' @return A data frame with metadata for each biopsy in GAMBL
#' @export
#'
#' @examples
#' #basic usage
#' my_metadata = get_gambl_metadata()
get_gambl_metadata = function(db="gambl_test",seq_type_filter = "genome",tissue_status_filter="tumour"){
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  sample_meta = tbl(con,"sample_metadata") %>% filter(seq_type == seq_type_filter & tissue_status == tissue_status_filter)
  #if we only care about genomes, we can drop/filter anything that isn't a tumour genome
  #The key for joining this table to the mutation information is to use sample_id. Think of this as equivalent to a library_id. It will differ depending on what assay was done to the sample.
  biopsy_meta = tbl(con,"biopsy_metadata") %>% select(-patient_id) %>% select(-pathology) %>% select(-time_point) %>% select(-EBV_status_inf) #drop duplicated columns
  all_meta = left_join(sample_meta,biopsy_meta,by="biopsy_id") %>% as.data.frame()

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
  DBI::dbDisconnect(con)
  return(all_meta)
}

#' Get the patient-centric clinical metadata
#'
#' @param db
#' @param time_unit Return follow-up times in one of three time units: year, month or day
#' @param censor_cbioportal Optionally request the censoring to be encoded in the specific style required by cBioPortal
#'
#' @return Data frame with one row for each patient_id
#' @export
#'
#' @examples
#' outcome_df = get_gambl_outcomes()
get_gambl_outcomes = function(db="gambl_test",time_unit="year",censor_cbioportal=FALSE){
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  all_outcome = tbl(con,"outcome_metadata") %>% as.data.frame()
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
  if(grepl("chr",chromosome)){
    chromosome = gsub("chr","",chromosome)
  }
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  if(!missing(region) || !missing(chromosome)){
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
#'
#' @examples
#' # basic usage
#' my_segments=get_cn_segments(region="chr8:128,723,128-128,774,067")
#' # specifying chromosome, start and end individually
#' my_segments=get_cn_segments(chromosome="8",qstart=128723128,qend=128774067)
#' # Asking for chromosome names to have a chr prefix (default is un-prefixed)
#' prefixed_segments = get_cn_segments(get_cn_segments(chromosome="12",qstart=122456912,qend=122464036,with_chr_prefix = TRUE))
get_cn_segments = function(db="gambl_test",table_name="seg_battenberg_hg19",chromosome=NULL,qstart=NULL,qend=NULL,region="",with_chr_prefix=FALSE){
  if(!region==""){
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
    filter(chrom == chromosome & start <= qstart & end >= qend) %>%
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
#'
#' @examples
#' #basic usage
#' my_mutations=get_ssm_by_region(region="chr8:128,723,128-128,774,067")
#' #specifying chromosome, start and end individually
#' my_mutations=get_ssm_by_region(chromosome="8",qstart=128723128,qend=128774067)
get_ssm_by_region = function(db="gambl_test",table_name="maf_slms3_hg19_icgc",chromosome,qstart,qend,region="",basic_columns=TRUE){
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
  if(basic_columns){
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
#'
#' @examples
#' #basic usage
#' maf_data = get_coding_ssm(limit_cohort=c("BL_Adult","BL_Pediatric","BL_ICGC"))
get_coding_ssm = function(db="gambl_test",table_name="maf_slms3_hg19_icgc",limit_cohort,exclude_cohort,limit_pathology,basic_columns=TRUE){
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Silent","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")
  sample_meta = tbl(con,"sample_metadata") %>% filter(seq_type == "genome" & tissue_status == "tumour")
  biopsy_meta = tbl(con,"biopsy_metadata") %>% select(-patient_id) %>% select(-pathology) %>% select(-time_point) %>% select(-EBV_status_inf) #drop duplicated columns
  all_meta = left_join(sample_meta,biopsy_meta,by="biopsy_id") %>% as.data.frame()

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
    filter(Variant_Classification %in% coding_class & Tumor_Sample_Barcode %in% sample_ids)

  muts = as.data.frame(muts)
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

