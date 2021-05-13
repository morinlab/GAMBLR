
maf_to_custom_track = function(maf_data){

}


#' Bring together all derived sample-level results from many GAMBL pipelines
#'
#' @return A table keyed on biopsy_id that contains a bunch of per-sample results from GAMBL
#' @export
#' @import tidyverse config
#'
#' @examples
collate_results = function(sample_table,write_to_file=FALSE){
  # important: if you are collating results from anything but WGS (e.g RNA-seq libraries) be sure to use biopsy ID as the key in your join
  # the sample_id should probably not even be in this file if we want this to be biopsy-centric
  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>% dplyr::select(sample_id,patient_id,biopsy_id)
  }
  #edit this function and add a new function to load any additional results into the main summary table

  sample_table = collate_sv_results(sample_table=sample_table)
  sample_table = collate_curated_sv_results(sample_table=sample_table)
  sample_table = collate_ashm_results(sample_table=sample_table)
  sample_table = collate_nfkbiz_results(sample_table=sample_table)
  sample_table = collate_sbs_results(sample_table=sample_table)
  sample_table = collate_derived_results(sample_table=sample_table)
  if(write_to_file){
    output_file = config::get("table_flatfiles")$derived
    output_base = config::get("repo_base")
    output_file = paste0(output_base,output_file)
    write_tsv(sample_table,file=output_file)
  }
  return(sample_table)
}

#' Extract derived results stored in the database (these are usually slower to derive on the fly)
#'
#' @param sample_table A data frame with sample_id as the first column
#'
#' @return Data frame with one row per sample. Contains the contents of the derived_data table in the database
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
collate_derived_results = function(sample_table){

  database_name = config::get("database_name")

  con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  derived_tbl = dplyr::tbl(con,"derived_data") %>% as.data.frame()
  derived_tbl = derived_tbl %>% dplyr::select(where( ~!all(is.na(.x)))) #drop the columns that are completely empty
  sample_table = dplyr::left_join(sample_table,derived_tbl)
  return(sample_table)
}

#' Title
#'
#' @param sample_table A data frame with sample_id as the first column
#' @param path_to_files Full local (base) path to the home of GAMBL outputs
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
collate_curated_sv_results = function(sample_table){
  path_to_files = config::get("derived_and_curated")
  project_base = config::get("project_base")
  #  "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
  manual_files = dir(paste0(project_base,path_to_files),pattern=".tsv")
  for(f in manual_files){
    full = paste0(project_base,path_to_files,f)
    this_data = read_tsv(full,comment = "#")
    #TO DO: fix this so it will join on biopsy_id or sample_id depending on which one is present
    sample_table = left_join(sample_table,this_data)
  }
  return(sample_table)
}

#' Annotate mutations with their copy number information
#'
#' @param this_sample Sample ID of the sample you want to annotate
#' @param coding_only Optional. set to TRUE to rescrict to only coding variants
#' @param from_flatfile Optional. instead of the database, load the data from a local MAF and seg file
#'
#' @return A list containing a data frame (MAF-like format) with two extra columns:
#' log.ratio is the log ratio from the seg file (NA when no overlap was found)
#' as well as the segmented copy number data with the same copy number information
#' CN is the rounded absolute copy number estimate of the region based on log.ratio (NA when no overlap was found)
#' @export
#' @import tidyverse data.table RMariaDB DBI dbplyr
#'
#' @examples
#' cn_list = assign_cn_to_ssm(this_sample="HTMCP-01-06-00422-01A-01D",coding_only=TRUE)
assign_cn_to_ssm = function(this_sample,coding_only=FALSE,from_flatfile=FALSE,
                            use_augmented_maf=FALSE){

  database_name = config::get("database_name")
  project_base = config::get("project_base")
  tool_name=config::get("analyses")$matched$copy_number


  #project_base = "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
  coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")
  if(from_flatfile){
    #get the genome_build for this sample
    bam_info = get_bams(this_sample)
    genome_build = bam_info$genome_build
    unix_group = bam_info$unix_group
    #maf path for a single file is easy to predict. This really should be generalized for all tools
    if(use_augmented_maf==TRUE){
      #results/gambl/rainstorm_circos/genome--grch37/01-augment_ssm/13-38657_tumorA--13-38657_normal--matched_slms-3.final_augmented.maf
      maf_path = paste0(project_base,unix_group,"/","rainstorm_circos/genome--",genome_build,"/01-augment_ssm/")
      this_sample_maf = dir(maf_path,pattern=paste0(this_sample,"--"))
      this_sample_maf = grep(".maf",this_sample_maf,value=T)
      this_sample_maf=paste0(maf_path,this_sample_maf)
    }else{
      slms3_path = paste0(project_base,unix_group,"/","slms-3_vcf2maf_current/99-outputs/genome--",genome_build,"/")
      this_sample_mafs = dir(slms3_path,pattern=paste0(this_sample,"--"))
      #use the lifted or native?
      this_sample_maf = this_sample_mafs[grep("converted",this_sample_mafs,invert=T)]
      this_sample_maf = paste0(slms3_path,this_sample_maf)
    }
    if(length(this_sample_maf)>1){
      print("WARNING: more than one MAF found for this sample. This shouldn't happen!")
      this_sample_maf = this_sample_maf[1]
    }
    #now we can load it
    maf_sample = fread_maf(this_sample_maf)

  }else{

    #get all the segments for a sample and filter the small ones then assign CN value from the segment to all SSMs in that region
    con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    maf_table = config::get("results_tables")$ssm
    maf_sample <- dplyr::tbl(con, maf_table) %>%
      dplyr::filter(Tumor_Sample_Barcode == this_sample) %>%
      as.data.frame()
  }
  if(coding_only){
    maf_sample = dplyr::filter(maf_sample,Variant_Classification %in% coding_class)
  }
  if(tool_name == "battenberg"){
    if(from_flatfile){
      battenberg_files = fetch_output_files(build=genome_build,base_path = "gambl/battenberg_current",tool="battenberg",search_pattern = ".igv.seg")
      battenberg_file = dplyr::filter(battenberg_files,tumour_sample_id==this_sample) %>%
        dplyr::pull(full_path) %>% as.character()
      if(length(battenberg_file)>1){
        print("WARNING: more than one SEG found for this sample. This shouldn't happen!")
        battenberg_file = battenberg_file[1]
      }
      seg_sample = read_tsv(battenberg_file) %>% as.data.table() %>% dplyr::mutate(size=end - start) %>%
        dplyr::filter(size > 100) %>%
        dplyr::mutate(chrom = gsub("chr","",chrom)) %>%
        dplyr::rename(Chromosome=chrom,Start_Position=start,End_Position=end)
    }else{
      seg_table = config::get("results_tables")$copy_number
      seg_sample = dplyr::tbl(con,seg_table) %>%
      dplyr::filter(ID == this_sample) %>%
      data.table::as.data.table() %>% dplyr::mutate(size=end - start) %>%
      dplyr::filter(size > 100) %>%
      dplyr::mutate(chrom = gsub("chr","",chrom)) %>%
      dplyr::rename(Chromosome=chrom,Start_Position=start,End_Position=end)
    }
    data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    a = data.table::as.data.table(maf_sample)
    a.seg = data.table::foverlaps(a, seg_sample, type="any")
    a$log.ratio = a.seg$log.ratio
    a$LOH = factor(a.seg$LOH_flag)
    a = dplyr::mutate(a,CN=round(2*2^log.ratio))
    seg_sample = dplyr::mutate(seg_sample,CN=round(2*2^log.ratio))
    seg_sample$LOH_flag = factor(seg_sample$LOH_flag)
    mutate(a,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
      #geom_point(aes(x=Start_Position,y=vaf,colour=CN),size=0.1)  +
      geom_segment(data=seg_sample,aes(x=Start_Position,xend=End_Position,y=CN,yend=CN,colour=LOH_flag)) +
      facet_wrap(~Chromosome,scales="free_x")
    return(list(maf=a,seg=seg_sample))
    if(!from_flatfile){
      DBI::dbDisconnect(con)
    }
  }else{
    print("ERROR: missing a required argument")
  }
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
refresh_full_table = function(table_name,connection,file_path){
  table_data = read_tsv(file_path)
  dbWriteTable(con,table_name,table_data,overwrite=TRUE)
  print(paste("POPULATING table:",table_name,"USING path:",file_path))
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
referesh_metadata_tables = function(){

  con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  all_metadata_info = sanity_check_metadata()
  tables = pull(all_metadata_info,table)
  files = pull(all_metadata_info,file)
  #lazily use a for loop for this
  for(i in c(1:length(files))){
    refresh_full_table(tables[i],con,files[i])
  }
}

sanity_check_metadata = function(){
  #e.g. for biopsy_metadata
  #biopsy_table = config::get("tables")$biopsies
  #metadata_files = config::get("table_flatfiles")$biopsies
  cfg = config::get("tables")
  database_name = config::get("database_name")
  metadata_tables = tibble(key=names(cfg),table=cfg) %>% unnest_auto("table")
  cfg = config::get("table_flatfiles")
  metadata_files = tibble(key=names(cfg),file=cfg) %>% unnest_auto("file")
  all_metadata_info = left_join(metadata_tables,metadata_files)
  base_path = config::get("repo_base")
  all_metadata_info = all_metadata_info %>% mutate(file=paste0(base_path,file))
  all_metadata_df = all_metadata_info %>% column_to_rownames(var = "key")
  #all samples with different seq_type and protocol must have a unique sample_id
  sample_df = read_tsv(all_metadata_df["samples","file"])
  tumour_samples = sample_df %>% select(patient_id,sample_id,biopsy_id,seq_type,protocol) %>%
    dplyr::filter(!is.na(biopsy_id))
  n_samp_bio = tumour_samples %>% count() %>% pull(n)
  #2876 unique samples
  #check for any multiplicity of sample_id
  n_samp = tumour_samples %>% select(-biopsy_id) %>% unique() %>% count() %>% pull(n)
  #should be the same number as above
  if(!n_samp == n_samp_bio){
    print("ERROR! some biopsies appear to have the same sample_id/protocol combination")
  }

  return(all_metadata_info)
}


collate_extra_metadata= function(sample_table,file_path){
  file_path = "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/experiments/2021-04-21-Trios-MiXCR/trios_relapse_timing.tsv"
  extra_df = read_tsv(file_path)
  sample_table = left_join(sample_table,extra_df,by=c("sample_id"="biopsy_id"))
}

#' Bring in the results from mutational signature analysis
#'
#' @param sample_table A data frame with sample_id as the first column
#' @param file_path
#'
#' @return
#' @export
#' @import tidyverse
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

#' Determine which cases have NFKBIZ UTR mutations
#'
#' @return
#' @export
#' @import tidyverse
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

#' Determine the hypermutation status of a few genes
#'
#' @return
#' @export
#' @import tidyverse
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

#' Determine and summarize which cases have specific oncogene SVs
#'
#' @param sample_table A data frame with sample_id as the first column
#' @param tool
#' @param oncogenes
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
collate_sv_results = function(sample_table,tool="manta",oncogenes=c("MYC","BCL2","BCL6","CCND1","IRF4")){
  if(tool=="manta"){
    all_svs = get_manta_sv()
  }
  annotated_svs = annotate_sv(all_svs) %>% dplyr::filter(!is.na(partner))
  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>% dplyr::select(sample_id,patient_id,biopsy_id)
  }
  multiout <- function(df, annotated, tool, oncogene_name) {
    some_fusions = dplyr::filter(annotated,gene==all_of(oncogene_name)) %>%
      group_by(tumour_sample_id) %>% arrange(partner) %>% dplyr::filter(row_number()==1)
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

#' Get some colour schemes for annotating figures
#'
#' @param classification (optionally request only colours for pathology, lymphgen or copy_number)
#'
#' @return A named vector of colour codes for lymphgen classes and pathology
#' @export
#' @import tidyverse
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
  copy_number_colours=c(
    "nLOH"="#E026D7",
    "8"="#380015",
    "7"="#380015",
    "6"="#380015",
    "5"="#67001F",
    "4"="#B2182B",
    "3"="#D6604D",
    "2"="#ede4c7",
    "1"="#92C5DE",
    "0"="#4393C3"
  )
  pathology_colours = c(
    "DLBCL"="#479450",
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
    "HGBL"="#B23F52",
    "SCBC"="#8c9c90",
    "UNSPECIFIED"="#cfba7c"
  )
  if(classification == "copy_number"){
    return(copy_number_colours)
  }
  if(classification == "pathology"){
    return(pathology_colours)
  }
  else{
    all_colours=c(lymphgen_colours,pathology_colours)
    return(all_colours)
  }
}

#' Get full paths for bam files for a sample or patient
#'
#' @param sample Either provide sample_id or patient_id
#' @param patient Either provide sample_id or patient_id
#'
#' @return A list that contains the genome_build and an igv-friendly build (igv_build), a list of bam file paths for tumour, normal and mrna data
#' @export
#' @import tidyverse
#'
#' @examples
#'
#' this_sv = filter(annotate_sv(get_manta_sv()),partner=="HIST1H2BK")
#' #arbitrarily grab one SV
#' bam_details = get_bams(sample=this_sv[1,"tumour_sample_id"])
get_bams = function(sample,patient){
  meta = get_gambl_metadata(tissue_status_filter = c("tumour","normal"))
  meta_mrna = get_gambl_metadata(seq_type_filter = "mrna")
  #get all samples for this patient
  if(missing(patient)){
    patient = meta %>% dplyr::filter(sample_id==sample) %>% dplyr::pull(patient_id)
  }
  meta_patient = meta %>% dplyr::filter(patient_id == patient)
  meta_mrna_patient = meta_mrna %>% dplyr::filter(patient_id == patient)
  build = dplyr::pull(meta_patient,genome_build) %>% head(1)
  if(build == "hs37d5"){
    igv_build = "hg19"
  }else{
    igv_build = build
  }
  tumour_genome_bams = dplyr::filter(meta_patient,seq_type == "genome" & tissue_status == "tumour") %>%
    dplyr::pull(data_path)
  bam_details = list(igv_build=igv_build, genome_build=build, tumour_bams=tumour_genome_bams)
  normal_genome_bams = dplyr::filter(meta_patient,seq_type == "genome" & tissue_status == "normal") %>%
    dplyr::pull(data_path)
  unix_group = dplyr::filter(meta_patient,seq_type == "genome" & tissue_status == "tumour") %>%
    dplyr::pull(unix_group) %>% unique()
  bam_details$unix_group = unix_group
  if(length(normal_genome_bams)){
    bam_details$normal_genome_bams = normal_genome_bams
  }
  if(length(normal_genome_bams)){
    bam_details$normal_genome_bams = normal_genome_bams
  }
  rnaseq_bams = dplyr::filter(meta_mrna_patient,seq_type == "mrna") %>% dplyr::pull(data_path)
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
#' @import tidyverse SRAdb
#'
#' @examples
#' #IMPORTANT: you must be running IGV on the host that is running R and you need to have it listening on a port
#' # The simplest scenario is to run this command on a terminal (if using a Mac), assuming you are using R on gphost10 and you have a ssh config that routes gp10 to that host
#' # ssh -X gp10
#' # then launch IGV (e.e. from a conda installation):
#' # conda activate igv; igv &
#' # this_sv = annotated_sv %>% filter(gene=="ETV6")
#' # tumour_bam = get_bams(sample=this_sv$tumour_sample_id)
#' # make_igv_snapshot(chrom=this_sv$chrom2, start=this_sv$start2, end=this_sv$end2, sample_id=this_sv$tumour_sample_id,out_path="~/IGV_snapshots/")
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

