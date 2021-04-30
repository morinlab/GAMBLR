


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
assign_cn_to_ssm = function(this_sample,coding_only=FALSE,from_flatfile=FALSE){

  database_name = config::get("database_name")
  project_base = config::get("project_base")
  tool_name=config::get("copy_number_tool")


  #project_base = "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
  coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")
  if(from_flatfile){
    #get the genome_build for this sample
    bam_info = get_bams(this_sample)
    genome_build = bam_info$genome_build
    unix_group = bam_info$unix_group
    #maf path for a single file is easy to predict. This really should be generalized for all tools
    slms3_path = paste0(project_base,unix_group,"/","slms-3_vcf2maf_current/99-outputs/genome--",genome_build,"/")
    this_sample_mafs = dir(slms3_path,pattern=paste0(this_sample,"--"))
    #use the lifted or native?
    this_sample_maf = this_sample_mafs[grep("converted",this_sample_mafs,invert=T)]
    this_sample_maf = paste0(slms3_path,this_sample_maf)
    if(length(this_sample_maf)>1){
      print("WARNING: more than one MAF found for this sample. This shouldn't happen!")
      this_sample_maf = this_sample_maf[1]
    }
    #now we can load it
    maf_sample = fread_maf(this_sample_maf)

  }else{
    #get all the segments for a sample and filter the small ones then assign CN value from the segment to all SSMs in that region
    con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    maf_sample <- dplyr::tbl(con, table_name) %>%
      dplyr::filter(Tumor_Sample_Barcode == this_sample) %>%
      as.data.frame()
  }
  if(coding_only){
    maf_sample = dplyr::filter(maf_sample,Variant_Classification %in% coding_class)
  }
  if(tool_name == "battenberg"){
    if(from_flatfile){
      battenberg_files = fetch_output_files(genome_build=genome_build,base_path = "gambl/battenberg_current",tool_name="battenberg",search_pattern = ".igv.seg")
      battenberg_file = filter(battenberg_files,tumour_sample_id==this_sample) %>% pull(full_path) %>% as.character()
      if(length(battenberg_file)>1){
        print("WARNING: more than one SEG found for this sample. This shouldn't happen!")
        battenberg_file = battenberg_file[1]
      }
      seg_sample = read_tsv(battenberg_file) %>% as.data.table() %>% dplyr::mutate(size=end - start) %>%
        dplyr::filter(size > 100) %>%
        dplyr::mutate(chrom = gsub("chr","",chrom)) %>%
        dplyr::rename(Chromosome=chrom,Start_Position=start,End_Position=end)
    }else{
      seg_sample = dplyr::tbl(con,"seg_battenberg_hg19") %>%
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
    #a$LOH = factor(a.seg$LOH_flag)
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
  con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  tables = pull(all_metadata_info,table)
  files = pull(all_metadata_info,file)
  #lazily use a for loop for this
  for(i in c(1:length(files))){
    refresh_full_table(tables[i],con,files[i])
  }
}

#' Title
#'
#' @param file_paths A vector of full file paths, e.g. the output of dir
#' @param tool_name The tool or pipeline that generated the files (should be the same for all)
#' @param output_type The file type to distinguish different output file types from the same pipeline (e.g. seg, maf, ploidy)
#' @param tool_version Optional: provide the version of the pipeline or tool
#' @param unix_group The unix group (should be the same for all)
#' @param sample_ids A vector of sample_id the same length and in the same order as the file paths
#'
#' @return
#' @export
#'
#' @examples
assemble_file_details = function(file_paths,tool_name,unix_group,sample_ids,output_type="ploidy",is_production="yes"){
  #the gambl_files table contains
  #sample_id, unix_group, tool_name, tool_version, seq_type, genome_build, is_production, output_type, is_lifted_over, pairing_status
  # is_production, output_type, file_path, file_timestamp

  files_df = tibble(sample_id = sample_ids,unix_group=unix_group,tool_name=tool_name,output_type=output_type,is_production=is_production,
                    file_timestamp = lapply(file_paths,function(x){as.character(file.mtime(x))}),
                    file_path=file_paths) %>% unnest_longer(file_timestamp)


  #date_info = file.mtime(file_path)
  con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  dbWriteTable(con,"gambl_files",files_df,append=TRUE)
}

#' Populate the database with the per-sample summarized results of various tools
#'
#' @param sample_table A data frame with sample_id as the first column
#' @param tool Name of the tool to get the results for
#' @param base_directory_gambl
#' @param base_directory_other
#'
#' @return Nothing
#' @export
#' @import tidyverse DBI
#'
#' @examples
populate_tool_results = function(sample_table){
  copy_number_tool = unlist(strsplit(config::get("copy_number_tool"),","))
  database_name = config::get("database_name")
  ssm_tool = config::get("ssm_tool")
  genome_build = unlist(strsplit(config::get("genome_builds"),","))
  for(gb in genome_buid){
    populate_each_tool_result(copy_number_tool,database_name,genome_build)
    populate_each_tool_result(ssm_tool,database_name,genome_build)
  }
}

#' Title
#'
#' @param tool
#' @param genome_build
#' @param database_name
#'
#' @return
#' @export
#'
#' @examples
populate_each_tool_result = function(tool,genome_build,database_name){
  con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  generic_update = function(field_name,sample_id,field_value){
    #note: we'll need to handle strings differently here once we start adding them
    for(i in c(1:length(field_value))){
      fv = field_value[i]
      sid = sample_id[i]
      if(is.numeric(fv)){
        update_q = paste0("UPDATE derived_data set ",field_name," = ", fv, " WHERE sample_id = \"", sid,"\";")
        print(update_q)
      }else{
        #need to add quotes to the value also
        update_q = paste0("UPDATE derived_data set ",field_name," = \"", fv, "\" WHERE sample_id = \"", sid,"\";")
        print(update_q)
      }
      dbExecute(con, update_q)
    }
  }
  #check if we're missing sample_ids from sample_table
  sample_table = config::get("tables")$samples
  derived_table = config::get("tables")$derived
  sample_ids = pull(sample_table,sample_id)
  for(id in sample_ids){
    check_q = paste0("select count(*) from ",derived_table, " where sample_id = \"",id,"\";")
    num = dbGetQuery(con,check_q)
    if(num ==0){
      insert_q = paste0("insert into ", derived_table, " (sample_id) values(\"",id,"\");")
      print(insert_q)
      dbExecute(con, insert_q)
    }
  }
  if(tool=="QC"){
    #flag any cases with QC issues raised
    collated = collate_results()
    qc_issues = select(collated,sample_id,QC_flag) %>% filter(!is.na(QC_flag)) %>% dplyr::mutate(flagged="yes")
    generic_update(sample_id=qc_issues$sample_id,field_name = "QC_issue",field_value = qc_issues$QC_flag)
  }
  if(tool == "sequenza"){
    parse_sequenza = function(sequenza_files){
      seq_data=sequenza_files %>%
        map(read_tsv) %>% #read each file into a list of tibbles
        map(head,1) %>% #just keep the first line
        reduce(rbind) %>% #rbind the elements all back into one
        rename(sequenza_cellularity=cellularity,sequenza_ploidy=ploidy) #change the column names
      return(seq_data)
    }

    seq_files_gambl_hg38 = fetch_output_files(genome_build="hg38",base_path = "gambl/sequenza_current",results_dir="02-sequenza",tool_name="sequenza")
    sequenza_results = parse_sequenza(seq_files_gambl_hg38$full_path)
    generic_update(sample_id=seq_files_gambl_hg38$tumour_sample_id,field_name="sequenza_purity",field_value=sequenza_results$sequenza_cellularity)
    generic_update(sample_id=seq_files_gambl_hg38$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)


    seq_files_gambl = fetch_output_files(genome_build="grch37",base_path = "gambl/sequenza_current",results_dir="02-sequenza",tool_name="sequenza")
    sequenza_results = parse_sequenza(seq_files_gambl$full_path)
    generic_update(sample_id=seq_files_gambl$tumour_sample_id,field_name="sequenza_purity",field_value=sequenza_results$sequenza_cellularity)
    generic_update(sample_id=seq_files_gambl$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)


    seq_files_gambl_hg38 = fetch_output_files(genome_build="hg38",base_path = "icgc_dart/sequenza_current",results_dir="02-sequenza",tool_name="sequenza")
    sequenza_results = parse_sequenza(seq_files_gambl_hg38$full_path)
    generic_update(sample_id=seq_files_gambl_hg38$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)

    seq_files_gambl_grch37 = fetch_output_files(genome_build="hs37d5",base_path = "icgc_dart/sequenza_current",results_dir="02-sequenza",tool_name="sequenza")
    sequenza_results = parse_sequenza(seq_files_gambl_grch37$full_path)
    generic_update(sample_id=seq_files_gambl_grch37$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)


  }
  if(tool == "slms3"){
    gambl_mut_maf <- tbl(con, "maf_slms3_hg19_icgc")
    #additional bookkeeping: set matched/unmatched information in the analysis table based on the matched normal ID
    gambl_mutation_normals = gambl_mut_maf %>% select(Tumor_Sample_Barcode) %>%
      group_by(Tumor_Sample_Barcode) %>% as.data.frame()
    gambl_meta_normals = get_gambl_metadata(tissue_status_filter=c('tumour','normal')) %>%
      select(patient_id,sample_id,tissue_status) %>%
      pivot_wider(id_cols=patient_id,names_from=tissue_status,values_from=sample_id)
    #the above has tumour and normal as separate columns for all paired samples.

    gambl_normals = mutate(gambl_normals,slms3_pairing_status= case_when(

    ))
    #just use the mutation table to get summary counts per sample and add to the derived table for convenience

    #update for gambl cases then do the same for icgc, then repeat for coding changes
    gambl_counts = gambl_mut_maf %>% select(Tumor_Sample_Barcode) %>%
      group_by(Tumor_Sample_Barcode) %>% tally() %>% as.data.frame()
    generic_update(sample_id=gambl_counts$Tumor_Sample_Barcode,field_name="slms3_ssm_total",field_value=gambl_counts$n)

    #gambl_vafs = gambl_maf %>% group_by(Tumor_Sample_Barcode) %>% mutate(vaf=mean(t_alt_count/(t_ref_count + t_alt_count))) %>%
    #  select(Tumor_Sample_Barcode, vaf) %>%
    #  as.data.frame()
    vaf_q = "select Tumor_Sample_Barcode, avg(t_alt_count/(t_ref_count + t_alt_count)) as vaf from maf_slms3_hg19_icgc group by Tumor_Sample_Barcode"

    vaf_tbl = dbGetQuery(con,vaf_q)

    #update all vafs at once
    generic_update(sample_id=vaf_tbl$Tumor_Sample_Barcode,field_name="slms3_mean_vaf",field_value=vaf_tbl$vaf)

    coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")

    #this is where things get SLOW! It seems the most efficient here to do this query per sample.
    for(sam in gambl_counts$Tumor_Sample_Barcode){
    #for(sam in some){
      #check here to see if there's a non-null value and skip if possible
      check_q = paste0("select count(*) as n from derived_data where sample_id = \"",sam,"\" and slms3_ssm_coding is not NULL;")
      num = dbGetQuery(con,check_q) %>% pull(n)
      print(paste(sam,num))
      if(num == 0){
        print(paste("working on:",sam))
        coding_num = gambl_mut_maf %>% filter(Tumor_Sample_Barcode == sam & Variant_Classification %in% coding_class) %>% count() %>% pull(n)
        generic_update(sample_id=sam,field_name="slms3_ssm_coding",field_value=coding_num)
      }else{
        print(paste0("skipping ", sam))
      }
    }
  }
  if(tool == "battenberg-ploidy"){
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
    "HGBL"="#B23F52"
  )
  if(classification == "copy_number"){
    print("copy number colours")
    return(copy_number_colours)
  }
  if(classification == "pathology"){
    print("pathology colours")
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
  unix_group = filter(meta_patient,seq_type == "genome" & tissue_status == "tumour") %>% pull(unix_group) %>% unique()
  bam_details$unix_group = unix_group
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

#' Use liftOver to convert a bedpe file between the two main genome builds (grch37/hg38)
#'
#' @param bedpe_file Either specify the path to a bedpe file
#' @param bedpe_df Or specify the bedpe data in a data frame
#' @param target_build Specify which build the data should be lifted to (must be one of hg19, grch37, hg38, grch38)
#'
#' @return Data frame containing original bedpe data with new coordinates
#' @export
#' @import tidyverse
#'
#' @examples
#' hg19_sv = get_manta_sv() %>% head(100)
#' hg38_sv = liftover_bedpe(bedpe_df=hg19_sv,target_build="hg38")
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
  char_vec = original_bedpe %>% tidyr::unite(united,sep="\t") %>% dplyr::pull(united)
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
#' @import tidyverse
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
#' @param base_path Either the full or relative path to where all the results directories are for the tool e.g. "gambl/sequenza_current"
#' @param results_dir
#' @param seq_type
#' @param genome_build
#' @param search_pattern
#'
#' @return A data frame with one row per file and sample IDs parsed from the file name along with other GAMBL wildcards
#' @export
#' @import tidyverse
#'
#' @examples
fetch_output_files = function(tool_name,base_path,results_dir="99-outputs",seq_type="genome",genome_build="hg38",search_pattern="cellularity_ploidy.txt"){
  if(!grepl("^/",base_path)){
    project_base = config::get("project_base")
    #project_base = "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
    base_path = paste0(project_base,base_path)
  }

  print(paste0("using results in: ",results_path))
  #path may contain either directories or files named after the sample pair
  dir_listing = dir(results_path,pattern="--")
  #start a data frame for tidy collation of details
  dir_df = tibble(short_path=dir_listing)  %>% mutate(sample = strsplit(short_path,"--"))
  unnested_df = dir_df %>% unnest_wider(sample,names_sep="_") %>% rename(tumour_sample_id=sample_1,normal_sample_id=sample_2)
  #find file with search_pattern per directory and handle any missing files. This is a bit slow.
  #unnested_df = unnested_df %>% head() %>% mutate(output_file = dir(paste0(results_path,short_path),pattern=search_pattern))
  #This still fails when a matching file isn't found. No clue why this doesn't work
  if(tool_name=="sequenza"){
    results_path = paste0(base_path,"/",results_dir,"/",seq_type,"--",genome_build,"/")
    print("THIS CAN BE SLOW!")
    unnested_df = unnested_df  %>% mutate(full_path=paste0(results_path,short_path,"/filtered/sequenza_alternative_solutions.txt"))
    named=pull(unnested_df,full_path)
    found_files=  tibble(filename=lapply(named,file.exists)) %>%
      unnest_longer(filename)
    new_df = cbind(unnested_df,found_files) %>% filter(found_files==TRUE)
    print(head(new_df))
  }else if(tool_name == "battenberg"){
    results_path = paste0(base_path,"/",results_dir,"/seg/",seq_type,"--",genome_build,"/")
    all_files = dir(results_path,pattern=search_pattern)
    #extract tumour and normal ID
    all_tumours = unlist(lapply(all_files,function(x){tumour=unlist(strsplit(x,"--"))[1]}))
    all_normals = unlist(lapply(all_files,function(x){tumour=unlist(strsplit(x,"--"))[2]}))

    all_files = unlist(lapply(all_files,function(x){paste0(results_path,x)}))
    new_df= data.frame(tumour_sample_id=all_tumours,normal_sample_id=all_normals,full_path=all_files)
    new_df = mutate(new_df,normal_sample_id = gsub(normal_sample_id,pattern = "_subclones.igv.seg",replacement = ""))
  }else if(tool_name == "battenberg_ploidy"){
    results_path = paste0(base_path,"/",results_dir,"/",seq_type,"--",genome_build,"/")
    print("THIS CAN BE SLOW!")
    unnested_df = unnested_df  %>% mutate(full_path=paste0(results_path,short_path))
    named=pull(unnested_df,full_path)
    found_files=  tibble(filename=lapply(named,function(x){dir(x,pattern=search_pattern)[1]})) %>%
      unnest_longer(filename)
  }
  return(new_df)
}



#' Title
#'
#' @param maf_file_path
#'
#' @return a data table containing MAF data from a MAF file
#' @export
#'
#' @examples
fread_maf = function(maf_file_path){
  maf_dt = data.table::fread(
    file = maf_file_path,
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
  return(maf_dt)
}



#' Title
#'
#' @param mafs
#' @param sample_id
#' @param genes
#' @param show_noncoding
#' @param detail
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
plot_multi_timepoint = function(mafs,sample_id,genes,show_noncoding=FALSE,detail){
  tp = c("A","B","C")
  title = paste(sample_id,detail,sep="\n")
  i = 1
  for (i in c(1:length(mafs))){
    maf_file = mafs[i]
    time_point=tp[i]
    print(paste(maf_file,time_point))
  }
  if(length(mafs)==2){
    A.maf = fread_maf(mafs[1])
    B.maf = fread_maf(mafs[2])

    A.maf = A.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=1) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    B.maf = B.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=2) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    all.maf=rbind(A.maf,B.maf)
    if(show_noncoding){
      coding.maf = subset_regions(all.maf,shm_regions)
    }
    else{
      coding.maf = filter(all.maf,!Variant_Classification %in% c("Silent","RNA","IGR","Intron","5'Flank","3'Flank","5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = filter(coding.maf,(Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }
    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = filter(coding.maf,coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = filter(coding.maf,coord %in% B.zero.coords)
    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords),"category"]="not-A"
    coding.maf[which(coord %in% B.zero.coords),"category"]="not-B"



    #actually this is changed in eitehr direction, not just gained
    just_gained_lg_all = filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" )
    #just_gained_lg = filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point !=2) %>%
    #  mutate(time_point = time_point +0.4)

    just_gained_lg = filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point == 2 & VAF > 0) %>%
      mutate(time_point = time_point +0.4)
    print(just_gained_lg)
    just_trunk = filter(coding.maf, Hugo_Symbol %in% genes & category == "trunk" & time_point ==1) %>%
      mutate(time_point = time_point-0.4)
    ggplot(coding.maf,aes(x=time_point,y=VAF,group=coord,colour=category)) +
      geom_point() + geom_line(alpha=0.5) +
      geom_text_repel(data=just_gained_lg,aes(label=Hugo_Symbol),size=4,segment.linetype=0) +
      geom_text_repel(data=just_trunk,aes(label=Hugo_Symbol),size=4,segment.linetype=0) +
      ggtitle(title) +
      theme_minimal()
    return(TRUE)
  }
  if(length(mafs)==3){
    A.maf = fread_maf(mafs[1])
    B.maf = fread_maf(mafs[2])
    C.maf = fread_maf(mafs[3])
    A.maf = A.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=1) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    B.maf = B.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=2) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    C.maf = C.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=3) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    all.maf=rbind(A.maf,B.maf,C.maf)
    if(show_noncoding){
      coding.maf = subset_regions(all.maf,shm_regions)
    }
    else{
      coding.maf = filter(all.maf,!Variant_Classification %in% c("Silent","RNA","IGR","Intron","5'Flank","3'Flank","5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = filter(coding.maf,(Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }

    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = filter(coding.maf,coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = filter(coding.maf,coord %in% B.zero.coords)

    C.rows = which(coding.maf$time_point==3)
    C.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==3 & VAF == 0), "coord"]))
    C.zero = filter(coding.maf,coord %in% C.zero.coords)

    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords),"category"]="not-A"
    coding.maf[which(coord %in% B.zero.coords),"category"]="not-B"
    coding.maf[which(coord %in% C.zero.coords),"category"]="not-C"

    just_gained_lg_all = filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" )
    just_gained_lg = filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" & time_point ==3) %>%
      mutate(time_point = time_point +0.4)

    just_trunk = filter(coding.maf, Hugo_Symbol %in% lg & category == "trunk" & time_point ==1) %>%
      mutate(time_point = time_point-0.4)
    ggplot(coding.maf,aes(x=time_point,y=VAF,group=coord,colour=category)) +
      geom_point() + geom_line(alpha=0.3) +
      geom_text_repel(data=just_gained_lg,aes(label=Hugo_Symbol),size=4,segment.linetype=0) +
      geom_text_repel(data=just_trunk,aes(label=Hugo_Symbol),size=4,segment.linetype=0) +
      ggtitle(title) +
      theme_minimal()
    return(TRUE)
  }
  return(FALSE)
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
get_gene_expression = function(database_name,hugo_symbols,tidy_expression_data,metadata,join_with="mrna"){
  if(missing(database_name)){
    database_name = config::get("database_name")
  }
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
      pivot_wider(names_from=Hugo_Symbol,values_from=expression)
    expression_wider = left_join(metadata,expression_wider,by=c("sample_id"="genome_sample_id"))
  }
  return(expression_wider)
}

#Not meant to be used routinely
#' Title
#'
#' @return
#' @export
#'
#' @examples
tidy_gene_expression = function(){
  #read in the full matrix
  ex_matrix_file=config::get("ex_matrix_file")
  ex_matrix_full = read_tsv(ex_matrix_file)

  ex_tidy = pivot_longer(ex_matrix_full,-Hugo_Symbol,names_to="sample_id",values_to="expression")
  ex_tidy_nfkbiz = filter(ex_tidy,Hugo_Symbol=="NFKBIZ")
  all_samples = pull(ex_tidy_nfkbiz,sample_id) %>% unique
  #retrieve the full list of sample_id for RNA-seq libraries that have data in this matrix

  #pull the full metadata for all RNA-seq samples in GAMBL
  #under the hood this is a join of the sample and biopsy tables but subset for RNA-seq
  rna_meta = get_gambl_metadata(seq_type_filter = "mrna")
  #subset to just the ones in the matrix and keep only the relevant rows
  rm_dupes = c("08-15460_tumorA","05-32150_tumorA") #toss two duplicated cases AND FFPE_Benchmarking
  #rna_meta[which(rna_meta$sample_id %in% rm_dupes,"cohort"]= "FFPE_Benchmarking"
  rna_meta_existing = rna_meta %>% filter(sample_id %in% all_samples) %>% select(sample_id,patient_id,biopsy_id,protocol)


  selected_libraries = rna_meta_existing %>%
  group_by(biopsy_id) %>%
    # Take the biopsy_id with the longest string length (e.g PolyA vs. Ribodepletion)
    slice_max(str_length(protocol), n = 1,with_ties=FALSE) %>%
    ungroup()


  #set the canonical library per biopsy_id by picking the ribominus for cases with more than one
  #duplicated = rna_meta_existing %>% group_by(biopsy_id) %>% tally() %>% filter(n>1) %>% select(biopsy_id)
  #singleton = rna_meta_existing %>% group_by(biopsy_id) %>% tally() %>% filter(n==1) %>% select(-n)

  #selected_duplicated = left_join(duplicated,rna_meta_existing,by="biopsy_id") %>%
  #  arrange(desc(protocol)) %>% arrange(biopsy_id) %>% group_by(biopsy_id) %>% filter(row_number()==1) %>% ungroup()
  #selected_singleton = left_join(singleton,rna_meta_existing,by="biopsy_id") %>% ungroup()
  #selected_all=rbind(selected_duplicated,selected_singleton)
  #put everything back together
  #at this point I have 1206 sample_ids
  ex_tidy = ex_tidy %>% rename(mrna_sample_id = sample_id)

  ex_tidy = ex_tidy %>% filter(mrna_sample_id %in% selected_libraries$sample_id)
  rna_meta = rna_meta %>% select(sample_id,biopsy_id)
  ex_tidy_final = left_join(ex_tidy,rna_meta,by=c("mrna_sample_id"="sample_id"))
  #this still has the mrna sample ID. Need to add the tumour_sample_id from this biopsy (where available)
  genome_meta = get_gambl_metadata() %>% select(biopsy_id,sample_id,patient_id,ffpe_or_frozen) %>% rename("genome_sample_id" = "sample_id")
  #this gets the metadata in the same format but restricted to genome samples
  #REMOVE ANNOYING DUPLICATE GENOMES!
  duplicated_cases = genome_meta %>% group_by(biopsy_id) %>% tally() %>% filter(n>1)

  selected_genomes = genome_meta %>%
    group_by(biopsy_id) %>%
    # Take the biopsy_id with the longest string length (e.g frozen)
    slice_max(str_length(ffpe_or_frozen), n = 1,with_ties=FALSE) %>%
    ungroup() %>% select(-patient_id,-ffpe_or_frozen)

  #join to genome metadata based on biopsy_id (should be the same for RNA-seq and tumour genomes)
  ex_tidy_genome = left_join(ex_tidy_final,selected_genomes,by="biopsy_id")

  ex_tidy_genome = select(ex_tidy_genome,Hugo_Symbol,mrna_sample_id,expression,biopsy_id,genome_sample_id)
  #write the data back out for use by others and loading into the database.
  tidy_expression_file = config::get("tidy_expression_file")
  write_tsv(ex_tidy_genome,file=tidy_expression_file)


}
