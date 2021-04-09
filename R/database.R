
require("dbplyr")
require("tidyverse")


#' Retrieve all copy number segments from the GAMBL database that overlap with a single genomic coordinate range
#'
#' @param db The GAMBL database name
#' @param table_name The table we are querying
#' @param chromosome The chromosome you are restricting to
#' @param start Start coordinate of the range you are restricting to
#' @param end End coordinate of the range you are restricting to
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately
#'
#' @return A data frame containing all the MAF data columns (one row per mutation)
#' @export
#'
#' @examples
#' #basic usage
#' my_segments=get_cn_segments(region="chr8:128,723,128-128,774,067")
#' #specifying chromosome, start and end individually
#' my_segments=get_cn_segments(chromosome="8",start=128723128,end=128774067)
get_cn_segments = function(db="gambl_test",table_name="seg_battenberg_hg19",chromosome=NULL,qstart=NULL,qend=NULL,region=""){
  if(!region==""){
    region = gsub(",","",region)
    #format is chr6:37060224-37151701
    split_chunks = unlist(strsplit(region,":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2],"-"))
    start=startend[1]
    end=startend[2]
  }
  #chr prefix the chromosome. This isn't yet standardized so it's just a workaround for now.
  if(grepl("chr",chromosome)){

  }else{
    chromosome = paste0("chr",chromosome)
  }
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  all_segs = tbl(con,table_name) %>%
    filter(chrom == chromosome & start <= qstart & end >= qend) %>%
    as.data.frame()
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
  return(as.data.frame(muts_gene))
}

#' Retrieve all SSMs from the GAMBL database within a single genomic coordinate range
#'
#' @param db The GAMBL database name
#' @param table_name The table we are querying
#' @param chromosome The chromosome you are restricting to
#' @param start Start coordinate of the range you are restricting to
#' @param end End coordinate of the range you are restricting to
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately
#'
#' @return A data frame containing all the MAF data columns (one row per mutation)
#' @export
#'
#' @examples
#' #basic usage
#' my_mutations=get_ssm_by_region(region="chr8:128,723,128-128,774,067")
#' #specifying chromosome, start and end individually
#' my_mutations=get_ssm_by_region(chromosome="8",start=128723128,end=128774067)
get_ssm_by_region = function(db="gambl_test",table_name="maf_slms3_hg19_icgc",chromosome=NULL,start=NULL,end=NULL,region=""){
  if(!region==""){
    region = gsub(",","",region)
    #format is chr6:37060224-37151701
    split_chunks = unlist(strsplit(region,":"))
    chromosome = split_chunks[1]
    chromosome = gsub("chr","",chromosome)
    startend = unlist(strsplit(split_chunks[2],"-"))
    start=startend[1]
    end=startend[2]
  }
  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = db)
  muts_region = tbl(con,table_name) %>%
    filter(Chromosome == chromosome & Start_Position > start & Start_Position < end)
  return(as.data.frame(muts_region))
}

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

  return(all_meta)
}
