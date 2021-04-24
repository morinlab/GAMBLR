require(tidyverse)
require(rtracklayer)
require(SRAdb)




#' Title
#'
#' @return A table keyed on biopsy_id that contains a bunch of per-sample results from GAMBL
#' @export
#'
#' @examples
collate_results = function(sample_table,write_to_file=FALSE){
  # important: if you are collating results from anything but WGS (e.g RNA-seq libraries) be sure to use biopsy ID as the key in your join
  # the sample_id should probably not even be in this file if we want this to be biopsy-centric
  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>% select(sample_id,patient_id,biopsy_id)
  }
  #edit this function and add a new function to load any additional results into the main summary table

  sample_table = collate_sv_results(sample_table=sample_table)
  sample_table = collate_curated_sv_results(sample_table=sample_table)
  sample_table = collate_ashm_results(sample_table=sample_table)
  sample_table = collate_nfkbiz_results(sample_table=sample_table)
  sample_table = collate_sbs_results(sample_table=sample_table)

  if(write_to_file){
    output_file = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/data/metadata/gambl_sample_results.tsv"
    write_tsv(sample_table,file=output_file)
  }
  return(sample_table)
}

#' Title
#'
#' @param sample_table
#' @param path_to_files
#'
#' @return
#' @export
#'
#' @examples
collate_curated_sv_results = function(sample_table,path_to_files="icgc_dart/derived_and_curated_metadata/"){
  project_base = "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
  manual_files = dir(paste0(project_base,path_to_files),pattern=".tsv")
  for(f in manual_files){
    full = paste0(project_base,"icgc_dart/derived_and_curated_metadata/",f)
    this_data = read_tsv(full)
    #TO DO: fix this so it will join on biopsy_id or sample_id depending on which one is present
    sample_table = left_join(sample_table,this_data)
  }
  return(sample_table)
}

#' Title
#'
#' @param sample_table
#' @param tool
#' @param base_directory_gambl
#' @param base_directory_other
#'
#' @return
#' @export
#'
#' @examples
populate_tool_results = function(database_name="gambl_test",sample_table,tool="battenberg",base_directory_gambl,base_directory_other){
  library(RMariaDB)
  library(DBI)
  con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  generic_update = function(field_name,sample_id,field_value){
    #note: we'll need to handle strings differently here once we start adding them
    for(i in c(1:length(field_value))){

      fv = field_value[i]
      sid = sample_id[i]
      update_q = paste0("UPDATE derived_data set ",field_name," = ", fv, " WHERE sample_id = \"", sid,"\";")
      print(update_q)
      dbExecute(con, update_q)
    }

  }
  #check if we're missing sample_ids from sample_table
  sample_ids = pull(sample_table,sample_id)
  for(id in sample_ids){
    check_q = paste0("select count(*) from derived_data where sample_id = \"",id,"\";")
    num = dbGetQuery(con,check_q)
    if(num ==0){
      insert_q = paste0("insert into derived_data (sample_id) values(\"",id,"\");")
      print(insert_q)
      dbExecute(con, insert_q)
    }
  }

  if(tool == "battenberg"){
    # parse purity and ploidy values from copy number caller and add to database
    parse_batt = function(batt_file){
      batt_data =  batt_file %>%
        map(read_tsv) %>%
        reduce(rbind) %>%
          rename(battenberg_cellularity=cellularity,battenberg_ploidy=ploidy,battenberg_psi=psi)
      return(batt_data)
    }
    files_gambl_hg38 = fetch_output_files(genome_build="hg38",base_path = "gambl/battenberg_current",results_dir="02-battenberg")
    #table column structure is as follows:
    # {tool}_{variable} e.g. battenberg_ploidy and battenberg_purity
    results_table = files_gambl_hg38 %>% mutate(parse_batt(full_path))
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_psi",field_value=results_table$battenberg_psi)
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_ploidy",field_value=results_table$battenberg_ploidy)
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_purity",field_value=results_table$battenberg_cellularity)

    files_gambl_grch37 = fetch_output_files(genome_build="grch37",base_path = "gambl/battenberg_current",results_dir="02-battenberg")
    results_table = files_gambl_grch37 %>% mutate(parse_batt(full_path))
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_psi",field_value=results_table$battenberg_psi)
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_ploidy",field_value=results_table$battenberg_ploidy)
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_purity",field_value=results_table$battenberg_cellularity)

    files_icgc_hg38 = fetch_output_files(genome_build="hg38",base_path = "icgc_dart/battenberg_current",results_dir="02-battenberg")
    results_table = files_icgc_hg38 %>% mutate(parse_batt(full_path))
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_psi",field_value=results_table$battenberg_psi)
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_ploidy",field_value=results_table$battenberg_ploidy)
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_purity",field_value=results_table$battenberg_cellularity)
    #note special genome build for icgc
    files_icgc_grch37 = fetch_output_files(genome_build="hs37d5",base_path = "icgc_dart/battenberg_current",results_dir="02-battenberg")
    results_table = files_icgc_grch37 %>% mutate(parse_batt(full_path))
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_psi",field_value=results_table$battenberg_psi)
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_ploidy",field_value=results_table$battenberg_ploidy)
    generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_purity",field_value=results_table$battenberg_cellularity)

  }
}
collate_extra_metadata= function(sample_table,file_path){
  file_path = "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/experiments/2021-04-21-Trios-MiXCR/trios_relapse_timing.tsv"
  extra_df = read_tsv(file_path)
  sample_table = left_join(sample_table,extra_df,by=c("sample_id"="biopsy_id"))
}

#' Title
#'
#' @param sample_table
#' @param file_path
#'
#' @return
#' @export
#'
#' @examples
collate_sbs_results = function(sample_table,file_path){
  if(missing(file_path)){
    file_path = "/projects/rmorin_scratch/prasath_scratch/gambl/sigprofiler/gambl_hg38/02-extract/slms3.gambl.icgc.hg38.matched.unmatched/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt"
  }
  signatures = read.csv(file_path,sep="\t",header=1,row.names=1)
  rs=rowSums(signatures)
  cn=colnames(signatures)
  new_sig = signatures
  for(col in cn){
    scaled_vals = signatures[,col] / rs
    new_sig[,col]=scaled_vals
  }
  sbs1 = signatures[,"SBS1"] / rs
  sbs9 = signatures[,"SBS9"] / rs
  sbs8 = signatures[,"SBS8"] / rs
  sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs9=sbs9,sbs8=sbs8)

  sample_table = left_join(sample_table,sbs)
  return(sample_table)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
collate_nfkbiz_results = function(sample_table){
  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>% select(sample_id,patient_id,biopsy_id)
  }
  this_region="chr3:101578214-101578365"
  nfkbiz_ssm = get_ssm_by_region(region=this_region) %>% pull(Tumor_Sample_Barcode) %>% unique
  nfkbiz_sv = get_manta_sv(region=this_region) %>% pull(tumour_sample_id) %>% unique
  nfkbiz = unique(c(nfkbiz_ssm,nfkbiz_sv))
  sample_table$NFKBIZ_UTR = "NEG"
  sample_table[sample_table$sample_id %in% nfkbiz,"NFKBIZ_UTR"]= "POS"
  return(sample_table)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
collate_ashm_results = function(sample_table){
  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>% select(sample_id,patient_id,biopsy_id)
  }
  #just annotate BCL2, MYC and CCND1 hypermutation
  regions_df = data.frame(name=c("CCND1","BCL2","MYC"),
        region=c("chr11:69455000-69459900","chr18:60983000-60989000","chr8:128747615-128751834"))
  region_mafs = lapply(regions_df$region,function(x){get_ssm_by_region(region=x,streamlined = TRUE)})
  tibbled_data = tibble(region_mafs, region_name = regions_df$name)
  unnested_df = tibbled_data %>% unnest_longer(region_mafs)
  unlisted_df = mutate(unnested_df,start=region_mafs$Start_Position,sample_id=region_mafs$Tumor_Sample_Barcode) %>%
    select(start,sample_id,region_name)
  tallied = unlisted_df %>% group_by(sample_id,region_name) %>%
    tally() %>%
    pivot_wider(values_from=n,names_from=region_name,values_fill=0,names_prefix="ashm_")

  sample_table = left_join(sample_table,tallied)
}

#' Title
#'
#' @param sample_table
#' @param tool
#' @param oncogenes
#'
#' @return
#' @export
#'
#' @examples
collate_sv_results = function(sample_table,tool="manta",oncogenes=c("MYC","BCL2","BCL6","CCND1","IRF4")){
  if(tool=="manta"){
    all_svs = get_manta_sv()
  }
  annotated_svs = annotate_sv(all_svs) %>% filter(!is.na(partner))
  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>% select(sample_id,patient_id,biopsy_id)
  }
  multiout <- function(df, annotated, tool, oncogene_name) {
    some_fusions = filter(annotated,gene==all_of(oncogene_name)) %>%
      group_by(tumour_sample_id) %>% arrange(partner) %>% filter(row_number()==1)
    df = mutate(df, "{tool}_{oncogene_name}_sv" := case_when(
      sample_id %in% some_fusions$tumour_sample_id ~ "POS",
      TRUE ~ "NEG"
    ))
    some_fusions = some_fusions %>% select(tumour_sample_id,partner)  %>% mutate("{tool}_{oncogene_name}_partner" := partner) %>% select(-partner)
    df = left_join(df,some_fusions,by=c("sample_id"="tumour_sample_id"))
    return(df)
  }
  out_table = sample_table

  for(oncogene in oncogenes){
    out_table = multiout(out_table,annotated_svs,"manta",oncogene)
  }
  return(out_table)
}

#' Title
#'
#' @param classification
#'
#' @return
#' @export
#'
#' @examples
get_gambl_colours = function(classification="lymphgen"){
  lymphgen_colours = c(
    "EZB" = "#F37A20",
    "ST2" = "#E73325",
    "BN2" = "#EB7EB1",
    "MCD" = "#3B5FAC",
    "N1" =  "#7F3293",
    "COMPOSITE" = "#7E8083",
    "Other" = "#55B55E"
  )
  pathology_colours = c(
    "B-ALL"="#C1C64B",
    "BL"="#926CAD",
    "FL"="#EA8368",
    "CLL"="#889BE5",
    "MCL"="#721F0F",
    "MM"="#CC9A42",
    "B-cell unclassified"="#B581C6",
    "COMFL"="#8BBC98",
    "PBL" = "#E058C0",
    "DLBCL-BL-like"="#34C7F4",
    "HGBL"="#B23F52"
  )
  all_colours=c(lymphgen_colours,pathology_colours)
  return(all_colours)
}

#' Get full paths for bam files for a sample or patient
#'
#' @param sample Either provide sample_id or patient_id
#' @param patient Either provide sample_id or patient_id
#'
#' @return A list that contains the genome_build and an igv-friendly build (igv_build), a list of bam file paths for tumour, normal and mrna data
#' @export
#'
#' @examples
#' this_sv = annotated_sv %>% filter(partner=="HIST1H2BK") %>% head(1)
#' #arbitrarily grab a SV
#' bam_details = get_bams(sample=this_sv$tumour_sample_id)
get_bams = function(sample,patient){
  meta = get_gambl_metadata(tissue_status_filter = c("tumour","normal"))
  meta_mrna = get_gambl_metadata(seq_type_filter = "mrna")
  #get all samples for this patient
  if(missing(patient)){
    patient = meta %>% filter(sample_id==sample) %>% pull(patient_id)
  }
  meta_patient = meta %>% filter(patient_id == patient)
  meta_mrna_patient = meta_mrna %>% filter(patient_id == patient)
  build = pull(meta_patient,genome_build) %>% head(1)
  if(build == "hs37d5"){
    igv_build = "hg19"
  }else{
    igv_build = build
  }
  tumour_genome_bams = filter(meta_patient,seq_type == "genome" & tissue_status == "tumour") %>% pull(data_path)
  bam_details = list(igv_build=igv_build, genome_build=build, tumour_bams=tumour_genome_bams)
  normal_genome_bams = filter(meta_patient,seq_type == "genome" & tissue_status == "normal") %>% pull(data_path)
  if(length(normal_genome_bams)){
    bam_details$normal_genome_bams = normal_genome_bams
  }
  if(length(normal_genome_bams)){
    bam_details$normal_genome_bams = normal_genome_bams
  }
  rnaseq_bams = filter(meta_mrna_patient,seq_type == "mrna") %>% pull(data_path)
  if(length(rnaseq_bams)){
    bam_details$rnaseq_bams = rnaseq_bams
  }
  return(bam_details)
}


#' Load bams and generate an IGV screenshot for one or more regions
#'
#' @param bams Character vector containing the full path to one or more bam files
#' @param genome_build String specifying the genome build for the bam files provided
#' @param region Optionally specify the region as a single string (e.g. "chr1:1234-1235")
#' @param region_bed Optionally specify regions in bed format (column order is assumed)
#' @param padding Optionally specify a positive value to broaden the region around the specified position
#' @param chrom Optionally specify the region by specifying the chromosome, start and end (see below)
#' @param start Optionally specify the region by specifying the start
#' @param end Optionally specify the region by specifying the end
#' @param sample_id Specify the sample_id or any other string you want embedded in the file name
#' @param out_path Specify the output directory where the snapshot will be written
#' @param igv_port Specify the port IGV is listening on
#'
#' @return
#' @export
#'
#' @examples
#' #IMPORTANT: you must be running IGV on the host that is running R and you need to have it listening on a port
#' # The simplest scenario is to run this command on a terminal (if using a Mac), assuming you are using R on gphost10 and you have a ssh config that routes gp10 to that host
#' # ssh -X gp10
#' # then launch IGV (e.e. from a conda installation):
#' # conda activate igv; igv &
#' this_sv = annotated_sv %>% filter(gene=="ETV6")
#' tumour_bam = get_bam_path(sample=this_sv$tumour_sample_id)
#' make_igv_snapshot(chrom=this_sv$chrom2, start=this_sv$start2, end=this_sv$end2, sample_id=this_sv$tumour_sample_id,out_path="~/IGV_snapshots/")
make_igv_snapshot = function(bams,genome_build,region,padding=200,chrom,start,end,sample_id,out_path="/tmp/",igv_port=60506){
  sock= IGVsocket(port = igv_port)
  IGVclear(sock)
  if(missing(region)){
    region=paste0(chrom,":",start-padding,"-",end+padding)
  }
  #region="chr6:90885409-90885410"
  IGVgenome(sock,genome=genome_build)
  IGVgoto(sock,region)
  for(bam_file in bams){
    IGVload(sock,bam_file)
  }
  filename = paste(sample_id,region,"snapshot.png",sep="_")
  IGVsnapshot(sock,fname=filename,dirname=out_path)
  return(paste0(out_path,filename))
}

#' Use liftOver to convert a bedpe file between the two main genome builds (grch37/hg38)
#'
#' @param bedpe_file Either specify the path to a bedpe file
#' @param bedpe_df Or specify the bedpe data in a data frame
#' @param target_build Specify which build the data should be lifted to (must be one of hg19, grch37, hg38, grch38)
#'
#' @return Data frame containing original bedpe data with new coordinates
#' @export
#'
#' @examples
#' hg38_sv = lifover_bedpe(bedpe_df=hg19_sv,target_build="hg38")
liftover_bedpe = function(bedpe_file,bedpe_df,target_build="grch37"){
  if(!missing(bedpe_file)){
    original_bedpe = read_tsv(bedpe_file,comment = "##",col_types="cddcddccccccccccccccccc")
  }
  if(!missing(bedpe_df)){
    original_bedpe = bedpe_df
  }
  if(!grepl("chr",original_bedpe$CHROM_A)){
    #add chr prefix
    original_bedpe = original_bedpe %>% mutate(CHROM_A = paste0("chr",CHROM_A)) %>% mutate(CHROM_B = paste0("chr",CHROM_B))
  }
  char_vec = original_bedpe %>% unite(united,sep="\t") %>% pull(united)
  bedpe_obj <- rtracklayer::import(text=char_vec,format="bedpe")
  this_patient = colnames(original_bedpe)[23]
  this_normal = colnames(original_bedpe)[22]
  if(target_build == "grch37" | target_build == "hg19"){
    chain = rtracklayer::import.chain(system.file("extdata","hg38ToHg19.over.chain",package="GAMBLR"))
  }else if(target_build == "grch38" | target_build == "hg38"){
    chain = rtracklayer::import.chain(system.file("extdata","hg19ToHg38.over.chain",package="GAMBLR"))
  }
  colnames(original_bedpe)[1]="CHROM_A"
  original_columns = colnames(original_bedpe)

  first_sv_lifted = rtracklayer::liftOver(bedpe_obj@first,chain)
  second_sv_lifted = rtracklayer::liftOver(bedpe_obj@second,chain)
  no_problem = !((elementNROWS(first_sv_lifted) != 1) | (elementNROWS(second_sv_lifted) != 1))
  first_ok = subset(first_sv_lifted,no_problem)
  second_ok = subset(second_sv_lifted,no_problem)
  first_ok_df = rtracklayer::export(first_ok,format="bed") %>% read_tsv(col_names = c("CHROM_A","START_A","END_A","name_A","score_A","STRAND_A")) %>% select(-score_A) %>% select(-name_A)
  second_ok_df = rtracklayer::export(second_ok,format="bed") %>% read_tsv(col_names = c("CHROM_B","START_B","END_B","name_B","score_B","STRAND_B")) %>% select(-score_B) %>% select(-name_B)
  ok_bedpe = original_bedpe[no_problem,]
  kept_cols = ok_bedpe %>% select(-c("CHROM_A","START_A","END_A","CHROM_B","START_B","END_B","STRAND_A","STRAND_B"))
  fully_lifted = bind_cols(first_ok_df,second_ok_df,kept_cols) %>% select(all_of(original_columns))
  return(fully_lifted)
}

#' This is a helper function that is not meant to be used routinely
#'
#' @param bedpe_paths
#' @param pattern
#' @param out_dir
#'
#' @return
#'
#' @examples
read_merge_manta_with_liftover = function(bedpe_paths=c(),pattern="--matched",out_dir){
  to_merge = list()
  for(thispath in bedpe_paths){

    #sample_paths = dir(thispath,pattern=paste0("--",pattern)) #skip unmatched cases for now
    sample_paths = dir(thispath,pattern=pattern) #DEBUGGING
    print(sample_paths)
    #sample_paths = head(sample_paths,15) #for debugging
    for(sp in sample_paths){
      full_path = paste0(thispath,sp,"/somaticSV.bedpe")
      print(paste("working on ", full_path))


      if(grepl("hg38",full_path)){
        print("using liftOver")
        svbed = read_and_liftover_bedpe(full_path) #load and convert to grch37 coordinates
      }else{
        svbed=read_tsv(full_path,comment = "##",col_types="cddcddccccccccccccccccc")
      }

      this_patient = colnames(svbed)[23]
      this_normal = colnames(svbed)[22]
      out_file = paste0(out_dir,"/",this_patient,"--",this_normal,"--hg38Togrch37_sv.tsv")
      print(paste("writing output to",out_file))
      infos = pull(svbed,this_patient)
      infos_n = pull(svbed,this_normal)
      colnames(svbed)[c(1:6)]=c("CHROM_A","START_A","END_A","CHROM_B","START_B","END_B")
      #all_vafs = get.sv.vaf(infos)
      #svbed$VAF = as.numeric(all_vafs)
      svbed$VAF_tumour = sapply(infos,function(x){as.numeric(tail(unlist(strsplit(x,":")),1))})
      svbed$DP_tumour = sapply(infos,function(x){as.numeric(tail(unlist(strsplit(x,":")),2)[1])})
      svbed$VAF_normal = sapply(infos_n,function(x){as.numeric(tail(unlist(strsplit(x,":")),1))})
      svbed$DP_normal = sapply(infos_n,function(x){as.numeric(tail(unlist(strsplit(x,":")),2)[1])})
      svbed$SOMATIC_SCORE = sapply(svbed$INFO_A,function(x){as.numeric(tail(unlist(strsplit(x,"=")),1))})
      #filter on PASS, score, VAF

      #svbed_filt = svbed %>% filter( SCORE > minScore & FILTER == "PASS") %>%
      #  dplyr::select(c(chrom1,start1,end1,chrom2,start2,end2))
      svbed$tumour_sample_id = this_patient
      svbed$normal_sample_id = this_normal
      if(grepl("--unmatched",sp)){
        svbed$pair_status = "unmatched"
      }else{
        svbed$pair_status = "matched"
      }
      print(head(svbed))
      svbed$NAME = "."

      svbed = svbed %>% select(CHROM_A,START_A,END_A,CHROM_B,START_B,END_B,NAME,SOMATIC_SCORE,STRAND_A,STRAND_B,TYPE,FILTER,VAF_tumour,VAF_normal,DP_tumour,DP_normal,tumour_sample_id,normal_sample_id,pair_status)
      #remove chr prefix from both chromosome names
      svbed = svbed %>% mutate(CHROM_A = gsub("chr","",CHROM_A)) %>% mutate(CHROM_B = gsub("chr","",CHROM_B))
      write_tsv(svbed,out_file,col_names=FALSE)
      #to_merge[[this_patient]] = svbed
    }
  }
}

#' Title
#'
#' @param tool_name
#' @param base_path
#' @param results_dir
#' @param seq_type
#' @param genome_build
#' @param search_pattern
#'
#' @return A data frame with one row per file and sample IDs parsed from the file name along with other GAMBL wildcards
#' @export
#'
#' @examples
fetch_output_files = function(tool_name,base_path,results_dir="99-outputs",seq_type="genome",genome_build="hg38",search_pattern="cellularity_ploidy.txt"){
  if(!grepl("^/",base_path)){
    project_base = "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
    base_path = paste0(project_base,base_path)
  }
  results_path = paste0(base_path,"/",results_dir,"/",seq_type,"--",genome_build,"/")
  print(paste0("using results in: ",results_path))
  print("THIS CAN BE SLOW!")
  #path may contain either directories or files named after the sample pair
  dir_listing = dir(results_path,pattern="--")
  #start a data frame for tidy collation of details
  dir_df = tibble(short_path=dir_listing)  %>% mutate(sample = strsplit(short_path,"--"))
  unnested_df = dir_df %>% unnest_wider(sample,names_sep="_") %>% rename(tumour_sample_id=sample_1,normal_sample_id=sample_2)
  #find file with search_pattern per directory and handle any missing files. This is a bit slow.
  #unnested_df = unnested_df %>% head() %>% mutate(output_file = dir(paste0(results_path,short_path),pattern=search_pattern))
  #This still fails when a matching file isn't found. No clue why this doesn't work
  unnested_df = unnested_df  %>% mutate(full_path=paste0(results_path,short_path))
  named=pull(unnested_df,full_path)

  found_files=  tibble(filename=lapply(named,function(x){dir(x,pattern=search_pattern)[1]})) %>%
    unnest_longer(filename)


  new_df = cbind(unnested_df,found_files) %>% filter(!is.na(filename))
  new_df = mutate(new_df,full_path=paste0(full_path,"/",filename))

  return(new_df)
}

