
#' Get the details including file paths for the anticipated outputs from a pipeline or tool
#'
#' @param Optionally provide a data frame with all file details
#' @param tool_name
#' @param unix_group
#' @param filename_end_pattern Optionally specify a pattern to search for the files among a longer set of files in the outputs
#' @param update_db Set to TRUE to overwrite any existing rows in the table for this tool/unix_group combination
#'
#' @return
#' @export
#'
#' @examples
find_expected_outputs = function(targ_df,tool_name,unix_group,filename_end_pattern,update_db=FALSE,target_path){
  repo_base =config::get("repo_base")
  if(missing(target_path)){
    target_path = paste0(repo_base,"targets/",tool_name,"--",unix_group)
  }

  if(tool_name == "manta"){
    if(missing(targ_df)){
      filename_end_pattern=".somaticSV.bedpe"
      targ_df = read_tsv(target_path,col_names = c("file")) %>%
        dplyr::filter(str_detect(file,pattern = filename_end_pattern))
      targ_df = mutate(targ_df,file_path=paste0(repo_base,file)) %>%
        separate(file,sep="/",into=c("results","unix_group","tool_version","outputs","type","seq_genome","detail","filename"))
      targ_df = separate(targ_df,seq_genome,sep="--",into=c("seq_type","genome_build")) %>%
        dplyr::select(-results,-type,-detail,-outputs) %>%
        separate(filename,sep="--",into=c("tumour_sample_id","normal_sample_id","pairing_status")) %>%
        mutate(pairing_status = str_remove(pairing_status,filename_end_pattern)) %>%
        separate(tool_version,sep="-",into=c("tool_name","tool_version"))
    }
    targ_df = targ_df %>%
      mutate(file_timestamp=file.info(file_path)$mtime)
    targ_df$output_type = "bedpe"
    print(targ_df)
    print(target_path)
  }else if(tool_name=="gridss"){
    filename_end_pattern=".gridss_somatic_filtered.bedpe"
    targ_df = read_tsv(target_path,col_names = c("file")) %>%
      dplyr::filter(str_detect(file,pattern = filename_end_pattern))
    targ_df = mutate(targ_df,file_path=paste0(repo_base,file)) %>%
      separate(file,sep="/",into=c("results","unix_group","tool_version","outputs","type","seq_genome","detail","filename"))
    targ_df = separate(targ_df,seq_genome,sep="--",into=c("seq_type","genome_build")) %>%
      dplyr::select(-results,-type,-detail,-outputs) %>%
      separate(filename,sep="--",into=c("tumour_sample_id","normal_sample_id","pairing_status")) %>%
      mutate(pairing_status = str_remove(pairing_status,filename_end_pattern)) %>%
      separate(tool_version,sep="-",into=c("tool_name","tool_version")) %>%
      mutate(file_timestamp=file.info(file_path)$mtime)
    targ_df$output_type = "bedpe"
  }
  if(update_db){
    database_name = config::get("database_name")

    con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    table_name = config::get("tables")$files
    message(paste("updating",table_name,"in",database_name))
    #clear all files for this tool/unix_group combination
    update_q = paste0("DELETE from ", table_name, " WHERE tool_name = \"", tool_name,"\" and unix_group = \"", unix_group, "\" ;")
    print(update_q)
    something = dbReadTable(conn=con,table_name)
    summarized = something %>% group_by(unix_group) %>% tally()
    #merge with incoming data
    print(summarized)
    #to_write = rbind(something,targ_df)
    dbExecute(con, update_q,immediate=TRUE)
    dbWriteTable(con,table_name,targ_df,append=TRUE)

    #dbExecute(con, update_q,immediate=TRUE)
    #dbWriteTable(con,table_name,expected_manta,append=TRUE)
    DBI::dbDisconnect(con)
  }
  return(targ_df)
}

#' Populate the database with the per-sample summarized results of various tools
#'
#' @param tool Name of the tool to get the results for
#'
#' @return Nothing
#' @export
#' @import tidyverse DBI
#'
#' @examples
populate_tool_results = function(){
  #IMPORTANT TODO: This function should only ever work with samples that exist in the metadata
  # Perhaps it should report any excluded outputs in case they need to be deleted from the main output directories
  matched_analyses = unlist(config::get("analyses")$matched)
  print(matched_analyses)
  database_name = config::get("database_name")

  genome_builds = unlist(strsplit(config::get("genome_builds"),","))
  groups= unlist(strsplit(config::get("unix_groups"),","))
  for(analysis_type in names(matched_analyses)){
    tool_name = matched_analyses[analysis_type]
    message(paste("populating results for",tool_name))
    populate_each_tool_result(tool=tool_name,genome_builds,groups)
  }
}




#' Title
#'
#' @param tool
#' @param genome_build A list of all genome builds to process
#' @param unix_group A list of all unix groups to process
#'
#' @return
#' @export
#'
#' @examples
populate_each_tool_result = function(tool,genome_builds,unix_groups){
  database_name = config::get("database_name")
  con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  all_meta = get_gambl_metadata()
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
  #sample_table = config::get("tables")$samples
  #derived_table = config::get("tables")$derived
  #sample_ids = pull(sample_table,sample_id)
  #for(id in sample_ids){
  #  check_q = paste0("select count(*) from ",derived_table, " where sample_id = \"",id,"\";")
  #  num = dbGetQuery(con,check_q)
  #  if(num ==0){
  #    insert_q = paste0("insert into ", derived_table, " (sample_id) values(\"",id,"\");")
  #    print(insert_q)
  #    dbExecute(con, insert_q)
  #  }
  #}
  if(tool=="QC"){
    #flag any cases with QC issues raised
    collated = collate_results()
    qc_issues = dplyr::select(collated,sample_id,QC_flag) %>%
      dplyr::filter(!is.na(QC_flag)) %>% dplyr::mutate(flagged="yes")
    generic_update(sample_id=qc_issues$sample_id,field_name = "QC_issue",field_value = qc_issues$QC_flag)
  }
  if(tool == "sequenza"){
    parse_sequenza = function(sequenza_files){
      seq_data=sequenza_files %>%
        purrr::map(read_tsv) %>% #read each file into a list of tibbles
        purrr::map(head,1) %>% #just keep the first line
        purrr::reduce(rbind) %>% #rbind the elements all back into one
        dplyr::rename(sequenza_cellularity=cellularity,sequenza_ploidy=ploidy) #change the column names
      return(seq_data)
    }

    #seq_files_gambl_hg38 = fetch_output_files(genome_build="hg38",base_path = "gambl/sequenza_current",results_dir="02-sequenza",tool_name="sequenza")
    #sequenza_results = parse_sequenza(seq_files_gambl_hg38$full_path)
    #generic_update(sample_id=seq_files_gambl_hg38$tumour_sample_id,field_name="sequenza_purity",field_value=sequenza_results$sequenza_cellularity)
    #generic_update(sample_id=seq_files_gambl_hg38$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)

    #seq_files_gambl = fetch_output_files(genome_build="grch37",base_path = "gambl/sequenza_current",results_dir="02-sequenza",tool_name="sequenza")
    #sequenza_results = parse_sequenza(seq_files_gambl$full_path)
    #generic_update(sample_id=seq_files_gambl$tumour_sample_id,field_name="sequenza_purity",field_value=sequenza_results$sequenza_cellularity)
    #generic_update(sample_id=seq_files_gambl$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)

    for(unix_group in unix_groups){
      for(genome_build in genome_builds){
        seq_files_gambl = fetch_output_files(build=genome_build,unix_group=unix_group,
                                         base_path = paste0(unix_group,"/sequenza_current"),
                                         results_dir="02-sequenza",tool="sequenza")
        sequenza_results = parse_sequenza(seq_files_gambl$full_path)

        generic_update(sample_id=seq_files_gambl$tumour_sample_id,field_name="sequenza_purity",field_value=sequenza_results$sequenza_cellularity)
        generic_update(sample_id=seq_files_gambl$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)
      }
    }


    #seq_files_gambl_hg38 = fetch_output_files(genome_build="hg38",base_path = "icgc_dart/sequenza_current",results_dir="02-sequenza",tool_name="sequenza")
    #sequenza_results = parse_sequenza(seq_files_gambl_hg38$full_path)
    #generic_update(sample_id=seq_files_gambl_hg38$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)

    #seq_files_gambl_grch37 = fetch_output_files(genome_build="hs37d5",base_path = "icgc_dart/sequenza_current",results_dir="02-sequenza",tool_name="sequenza")
    #sequenza_results = parse_sequenza(seq_files_gambl_grch37$full_path)
    #generic_update(sample_id=seq_files_gambl_grch37$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)


  }
  if(tool == "manta"){
    #fetch the output file names per group/build combination
    message("processing results from manta")
    message(paste(unix_groups,sep=","))
    #separately process by unix group
    for(ug in unix_groups){
      #files_df = find_files_extract_wildcards(tool_name="manta",genome_build=genome_builds,search_pattern=".bed",unix_group=ug)
      files_df = find_expected_outputs(tool_name="manta",unix_group=ug)
      print(head(files_df))
      message(paste("processing",ug))
      n_missing =  files_df %>% dplyr::filter(is.na(file_timestamp)) %>% count() %>% pull(n)
      if(n_missing){
        message(paste("missing outputs for",n_missing))
        files_df=files_df %>% dplyr::filter(!is.na(file_timestamp))
      }
      #are there unexpected outputs?
      dupes = names(which(table(files_df$tumour_sample_id)>1))
      if(length(dupes)>0){
        message("DUPLICATE RESULTS EXIST. PLEASE FIX THIS AND RERUN")
        message(dupes)
        return()
      }

      all_tsb = pull(all_meta,sample_id)
      no_meta = unique(pull(files_df[which(!files_df$tumour_sample_id %in% all_tsb),],tumour_sample_id))
      if(length(no_meta)>0){
        message("DROPPING RESULTS FROM SAMPLES WITH NO METADATA:")
        print(no_meta)
        files_df = files_df %>% dplyr::filter(!tumour_sample_id %in% no_meta)
      }
      #TODO: Flag and drop files with too many SVs before merging (i.e. remove really bad data). Done manually currently.
      manta_df = process_all_manta_bedpe(files_df,group=ug) #need to add this to the database. Not currently automated
      manta_df %>% group_by(tumour_sample_id) %>% tally() %>% arrange(desc(n)) #have a look at the top offenders (most SVs)


    }
  }
  if(tool == "slms3"){
    gambl_mut_maf <- tbl(con, "maf_slms3_hg19_icgc")
    #additional bookkeeping: set matched/unmatched information in the analysis table based on the matched normal ID
    gambl_mutation_normals = gambl_mut_maf %>% dplyr::select(Tumor_Sample_Barcode) %>%
      group_by(Tumor_Sample_Barcode) %>% as.data.frame()
    gambl_meta_normals = get_gambl_metadata(tissue_status_filter=c('tumour','normal')) %>%
      dplyr::select(patient_id,sample_id,tissue_status) %>%
      pivot_wider(id_cols=patient_id,names_from=tissue_status,values_from=sample_id)
    #the above has tumour and normal as separate columns for all paired samples.

    gambl_normals = mutate(gambl_normals,slms3_pairing_status= case_when(

    ))
    #just use the mutation table to get summary counts per sample and add to the derived table for convenience

    #update for gambl cases then do the same for icgc, then repeat for coding changes
    gambl_counts = gambl_mut_maf %>% dplyr::select(Tumor_Sample_Barcode) %>%
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
        coding_num = gambl_mut_maf %>%
          dplyr::filter(Tumor_Sample_Barcode == sam & Variant_Classification %in% coding_class) %>%
          dplyr::count() %>% dplyr::pull(n)
        generic_update(sample_id=sam,field_name="slms3_ssm_coding",field_value=coding_num)
      }else{
        print(paste0("skipping ", sam))
      }
    }
  }
  if(tool == "battenberg_ploidy"){
    # parse purity and ploidy values from copy number caller and add to database
    parse_batt = function(batt_file){
      batt_data =  batt_file %>%
        purrr::map(read_tsv) %>%
        purrr::reduce(rbind) %>%
        dplyr::rename(battenberg_cellularity=cellularity,battenberg_ploidy=ploidy,battenberg_psi=psi)
      return(batt_data)
    }


    for(unix_group in unix_groups){
      message(unix_group)
      for(genome_build in genome_builds){
        message(genome_build)
        files = fetch_output_files(build=genome_build,unix_group=unix_group,
                                             base_path = paste0(unix_group,"/battenberg_current"),
                                             results_dir="02-battenberg",tool="battenberg_ploidy")

        results_table = files %>% mutate(parse_batt(full_path))
        generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_psi",field_value=results_table$battenberg_psi)
        generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_ploidy",field_value=results_table$battenberg_ploidy)
        generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_purity",field_value=results_table$battenberg_cellularity)
        #results_table = files %>% dplyr::mutate(parse_batt(full_path))
        #sequenza_results = parse_sequenza(seq_files_gambl$full_path)

        #generic_update(sample_id=seq_files_gambl$tumour_sample_id,field_name="sequenza_purity",field_value=sequenza_results$sequenza_cellularity)
        #generic_update(sample_id=seq_files_gambl$tumour_sample_id,field_name="sequenza_ploidy",field_value=sequenza_results$sequenza_ploidy)
      }
    }


    #files_gambl_hg38 = fetch_output_files(genome_build="hg38",base_path = "gambl/battenberg_current",results_dir="02-battenberg")
    #table column structure is as follows:
    # {tool}_{variable} e.g. battenberg_ploidy and battenberg_purity



    #files_gambl_grch37 = fetch_output_files(genome_build="grch37",base_path = "gambl/battenberg_current",results_dir="02-battenberg")
    #results_table = files_gambl_grch37 %>% mutate(parse_batt(full_path))
    #generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_psi",field_value=results_table$battenberg_psi)
    #generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_ploidy",field_value=results_table$battenberg_ploidy)
    #generic_update(sample_id=results_table$tumour_sample_id,field_name="battenberg_purity",field_value=results_table$battenberg_cellularity)


  }
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
  print(head(bedpe_paths))
  for(thispath in bedpe_paths){
    sample_paths = dir(thispath,pattern=pattern) #DEBUGGING
    print(sample_paths)
    #sample_paths = head(sample_paths,15) #for debugging
    for(sp in sample_paths){
      full_path = paste0(thispath,sp,"/somaticSV.bedpe")
      print(paste("working on HERE:", full_path))
      if(grepl("hg38",full_path)){
        print("using liftOver")
        svbed = liftover_bedpe(full_path) #load and convert to grch37 coordinates
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

      svbed = svbed %>% dplyr::select(CHROM_A,START_A,END_A,CHROM_B,START_B,END_B,NAME,SOMATIC_SCORE,STRAND_A,STRAND_B,TYPE,FILTER,VAF_tumour,VAF_normal,DP_tumour,DP_normal,tumour_sample_id,normal_sample_id,pair_status)
      #remove chr prefix from both chromosome names
      svbed = svbed %>% mutate(CHROM_A = gsub("chr","",CHROM_A)) %>% mutate(CHROM_B = gsub("chr","",CHROM_B))
      write_tsv(svbed,out_file,col_names=FALSE)
      #to_merge[[this_patient]] = svbed
    }
  }
}


#' This is a helper function that is not meant to be used routinely
#'
#' @param bedpe_paths
#' @param pattern
#' @param out_dir
#' @param projection_build The genome we want all results to be relative to (lifted if necessary)
#'
#' @return
#' @import tidyverse
#'
#' @export
#'
#' @examples
process_all_manta_bedpe = function(file_df,out_dir,group,genome_build,projection_build="grch37"){
  to_merge = list()
  if(missing(out_dir)){
    project_base =config::get("project_base")
    base_out_dir = config::get("results_staging")$manta
    out_dir = paste0(project_base,group,"/",base_out_dir)
  }

  process_manta = function(bedpe_file,liftover_to_hg19=FALSE,liftover_to_hg38=FALSE,only_return_missing=FALSE,projection="grch37"){
    cnames = c("CHROM_A","START_A","END_A","CHROM_B","START_B","END_B","NAME","SOMATIC_SCORE","STRAND_A","STRAND_B","TYPE","FILTER","VAF_tumour","VAF_normal","DP_tumour","DP_normal","tumour_sample_id","normal_sample_id","pair_status")

    svbed=read_tsv(bedpe_file,comment = "##",col_types="cddcddccccccccccccccccc")
    this_patient = colnames(svbed)[23]
    this_normal = colnames(svbed)[22]

    if(grepl("--unmatched",bedpe_file)){
      pair_status = "unmatched"
      svbed$pair_status = "unmatched"
    }else{
      svbed$pair_status = "matched"
      pair_status = "matched"
    }
    is_lifted = "native"
    if(liftover_to_hg19 || liftover_to_hg38){
      is_lifted = "lifted"
    }

    if(genome_build == projection | (genome_build == "hs37d5" & projection == "grch37")){
      is_lifted = "native"
    }
    out_file = paste0(out_dir,"/",this_patient,"--",this_normal,"--",pair_status,"--", is_lifted,"--genome--",genome_build,"--",projection_build,"_sv.tsv")
    message("working on OVER HERE:",bedpe_file)
    print(paste("output:",out_file))
    if(file.exists(out_file)){
      if(!only_return_missing){
        print(paste("LOADING",out_file))
        svbed = read_tsv(out_file,col_types = "ccccccccccccnnnnccc",col_names = cnames)
        return(svbed)
      }
      else{
        svbed = dplyr::filter(svbed,is.na(tumour_sample_id))
        return(svbed)
      }
    }
    if(liftover_to_hg19){
      svbed = liftover_bedpe(bedpe_df=svbed)
    }else if(liftover_to_hg38){
      svbed = liftover_bedpe(bedpe_df=svbed,target_build = "hg38")
    }

    infos = pull(svbed,this_patient)
    infos_n = pull(svbed,this_normal)
    colnames(svbed)[c(1:6)]=c("CHROM_A","START_A","END_A","CHROM_B","START_B","END_B")

    svbed$VAF_tumour = sapply(infos,function(x){as.numeric(tail(unlist(strsplit(x,":")),1))})
    svbed$DP_tumour = sapply(infos,function(x){as.numeric(tail(unlist(strsplit(x,":")),2)[1])})
    svbed$VAF_normal = sapply(infos_n,function(x){as.numeric(tail(unlist(strsplit(x,":")),1))})
    svbed$DP_normal = sapply(infos_n,function(x){as.numeric(tail(unlist(strsplit(x,":")),2)[1])})
    svbed$SOMATIC_SCORE = sapply(svbed$INFO_A,function(x){as.numeric(tail(unlist(strsplit(x,"=")),1))})

    svbed$tumour_sample_id = this_patient
    svbed$normal_sample_id = this_normal
    message(paste("checking status:",bedpe_file))


    svbed$NAME = "."
    svbed = svbed %>% dplyr::select(CHROM_A,START_A,END_A,CHROM_B,START_B,END_B,NAME,SOMATIC_SCORE,STRAND_A,STRAND_B,TYPE,FILTER,VAF_tumour,VAF_normal,DP_tumour,DP_normal,tumour_sample_id,normal_sample_id,pair_status)
    #remove chr prefix from both chromosome names
    svbed = svbed %>% mutate(CHROM_A = gsub("chr","",CHROM_A)) %>% mutate(CHROM_B = gsub("chr","",CHROM_B))
    #print(paste("writing output to",out_file))

    #run liftover after formatting?

    write_tsv(svbed,out_file,col_names=FALSE)
    return(svbed)
  }

  #separately run the hg38 and other builds, separately run per unix_group
  if(projection_build=="grch37"){
    if(genome_build == "hg38"){
      hg38_files = dplyr::filter(file_df,genome_build == "hg38" & unix_group == group) %>% pull(file_path)
      bed_data_lifted = hg38_files %>%
        purrr::map(process_manta,liftover_to_hg19=TRUE,projection=projection_build) %>%
        purrr::reduce(rbind)
    }else{
      not_hg38_files = dplyr::filter(file_df,genome_build != "hg38" & unix_group == group) %>% pull(file_path)
      bed_data_not_lifted =not_hg38_files %>%
        purrr::map(process_manta,liftover_to_hg19=FALSE,projection=projection_build) %>%
        purrr::reduce(rbind)
    }
  }else if(projection_build == "hg38"){
    if(genome_build == "hg38"){
      hg38_files = dplyr::filter(file_df,genome_build == "hg38" & unix_group == group) %>% pull(file_path)
      bed_data_not_lifted = hg38_files %>%
        purrr::map(process_manta,liftover_to_hg38=FALSE,liftover_to_hg19=FALSE,projection=projection_build) %>%
        purrr::reduce(rbind)
    }else{
      not_hg38_files = dplyr::filter(file_df,genome_build != "hg38" & unix_group == group) %>% pull(file_path)
      bed_data_lifted =not_hg38_files %>%
        purrr::map(process_manta,liftover_to_hg38=TRUE,liftover_to_hg19=FALSE,projection=projection_build) %>%
        purrr::reduce(rbind)
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
fetch_output_files = function(tool,unix_group,base_path,results_dir="99-outputs",
                              seq_type="genome",build="hg38",
                              search_pattern="cellularity_ploidy.txt"){
  if(!grepl("^/",base_path)){
    project_base = config::get("project_base")
    #project_base = "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
    base_path = paste0(project_base,base_path)
  }
  if(tool == "battenberg"){
    results_path = paste0(base_path,"/",results_dir,"/seg/",seq_type,"--",build,"/")
  }else{
    results_path = paste0(base_path,"/",results_dir,"/",seq_type,"--",build,"/")
  }
  #path may contain either directories or files named after the sample pair

  dir_listing = dir(results_path,pattern="--")
  #start a data frame for tidy collation of details
  dir_df = tibble(short_path=dir_listing)  %>% mutate(sample = strsplit(short_path,"--"))
  unnested_df = dir_df %>% unnest_wider(sample,names_sep="_") %>%
    dplyr::rename(tumour_sample_id=sample_1,normal_sample_id=sample_2)
  #find file with search_pattern per directory and handle any missing files. This is a bit slow.
  #unnested_df = unnested_df %>% head() %>% mutate(output_file = dir(paste0(results_path,short_path),pattern=search_pattern))
  #This still fails when a matching file isn't found. No clue why this doesn't work
  if(tool=="sequenza"){

    print(paste0("using results in: ",results_path))
    print("THIS CAN BE SLOW!")
    unnested_df = unnested_df  %>% mutate(full_path=paste0(results_path,short_path,"/filtered/sequenza_alternative_solutions.txt"))
    named=pull(unnested_df,full_path)
    found_files=  tibble(filename=lapply(named,file.exists)) %>%
      unnest_longer(filename)
    new_df = cbind(unnested_df,found_files) %>% dplyr::filter(found_files==TRUE)
    print(head(new_df))
  }else if(tool == "battenberg"){
    results_path = paste0(base_path,"/",results_dir,"/seg/",seq_type,"--",build,"/")
    all_files = dir(results_path,pattern=search_pattern)
    #extract tumour and normal ID
    all_tumours = unlist(lapply(all_files,function(x){tumour=unlist(strsplit(x,"--"))[1]}))
    all_normals = unlist(lapply(all_files,function(x){tumour=unlist(strsplit(x,"--"))[2]}))

    all_files = unlist(lapply(all_files,function(x){paste0(results_path,x)}))
    new_df= data.frame(tumour_sample_id=all_tumours,normal_sample_id=all_normals,full_path=all_files)
    new_df = mutate(new_df,normal_sample_id = gsub(normal_sample_id,pattern = "_subclones.igv.seg",replacement = ""))
  }else if(tool == "battenberg_ploidy"){
    #results_path = paste0(base_path,"/",results_dir,"/",seq_type,"--",build,"/")
    search_pattern="_cellularity_ploidy.txt"
    print("THIS CAN BE SLOW!")
    unnested_df = unnested_df  %>% mutate(full_path=paste0(results_path,short_path))
    named=pull(unnested_df,full_path)
    found_files=  tibble(filename=lapply(named,function(x){dir(x,pattern=search_pattern)[1]})) %>%
      unnest_longer(filename) #%>% mutate(path=paste0(full_path,filename))
    found_files = cbind(found_files,unnested_df) %>% dplyr::filter(!is.na(filename)) %>% mutate(full_path=paste0(full_path,"/",filename))
    return(found_files)
  }else if(tool == "manta"){
    tool_results_path = config::get("results_directories")$manta

    search_pattern=".bedpe"
    new_df = find_files_extract_wildcards(tool_name="manta",genome_build=c("hg38","grch37"),search_pattern=".bed")

    #details to include
    #assemble_file_details = function(file_paths,tool_name,unix_group,sample_ids,output_type="ploidy",is_production="yes")
  }
  return(new_df)
}



#' Title
#'
#' @param tool_results_path
#' @param search_pattern
#' @param genome_build
#' @param seq_type
#' @param unix_group
#'
#' @return
#' @export
#'
#' @examples
#' file_details_manta = find_files_extract_wildcards(tool_name="manta",genome_build=c("hg38","grch37"),search_pattern=".bed")
find_files_extract_wildcards = function(tool_results_path,search_pattern,genome_build,seq_type="genome",unix_group="gambl",tool_name){
  project_base = config::get("project_base")
  if(missing(tool_results_path)){
    tool_results_paths = config::get("results_directories")
    tool_results_path = tool_results_paths[[tool_name]]
  }
  results_paths = paste0(project_base,unix_group,"/",tool_results_path,"genome--",genome_build,"/somaticSV/")
  found_files=  tibble(filename=lapply(results_paths,function(x){dir(x,pattern=search_pattern)})) %>%
    mutate(path=results_paths,genome_build=genome_build) %>%
    unnest_longer(filename) %>% dplyr::filter(!is.na(filename)) %>%
    mutate(file_path=paste0(path,filename)) %>%
    mutate(pairing_status = case_when(
      grepl("--unmatched",filename) ~ "unmatched",
      TRUE ~"matched"
    )) %>%
    dplyr::mutate(sample_id = strsplit(filename,"--")) %>%
    unnest_wider(sample_id,names_sep = "-") %>% dplyr::rename(tumour_sample_id=`sample_id-1`,normal_sample_id=`sample_id-2`) %>%
    dplyr::mutate(tool_name=tool_name,seq_type=seq_type,unix_group=unix_group) %>%
    dplyr::select(tumour_sample_id,unix_group,tool_name,seq_type,genome_build,file_path,pairing_status,normal_sample_id)
  return(found_files)
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

#' Read a full expression matrix and subset to samples in GAMBL that have metadata (remove duplicates with consistent preferences)
#'
#' @return
#' @export
#'
#' @examples
tidy_gene_expression = function(){
  #read in the full matrix
  ex_matrix_file = config::get("results_merged")$ex_matrix_file
  tidy_expression_file = config::get("results_merged")$tidy_expression_file
  ex_matrix_full = read_tsv(ex_matrix_file)
  ex_matrix_full = dplyr::select(ex_matrix_full,-gene_id,-ensembl_gene_id) %>%
    dplyr::rename("Hugo_Symbol"="hgnc_symbol")

  ex_tidy = pivot_longer(ex_matrix_full,-Hugo_Symbol,names_to="sample_id",values_to="expression")
  ex_tidy_nfkbiz = dplyr::filter(ex_tidy,Hugo_Symbol=="NFKBIZ")
  all_samples = pull(ex_tidy_nfkbiz,sample_id) %>% unique
  #retrieve the full list of sample_id for RNA-seq libraries that have data in this matrix

  #pull the full metadata for all RNA-seq samples in GAMBL
  #under the hood this is a join of the sample and biopsy tables but subset for RNA-seq
  rna_meta = get_gambl_metadata(seq_type_filter = "mrna")
  #subset to just the ones in the matrix and keep only the relevant rows
  rm_dupes = c("08-15460_tumorA","05-32150_tumorA") #toss two duplicated cases AND FFPE_Benchmarking
  #rna_meta[which(rna_meta$sample_id %in% rm_dupes,"cohort"]= "FFPE_Benchmarking"
  rna_meta_existing = rna_meta %>%
    dplyr::filter(sample_id %in% all_samples) %>%
    dplyr::select(sample_id,patient_id,biopsy_id,protocol)


  selected_libraries = rna_meta_existing %>%
    group_by(biopsy_id) %>%
    # Take the biopsy_id with the longest string length (e.g PolyA vs. Ribodepletion)
    slice_max(str_length(protocol), n = 1,with_ties=FALSE) %>%
    ungroup()
  #this brings the total from 1399 to 1298

  #set the canonical library per biopsy_id by picking the ribominus for cases with more than one
  #duplicated = rna_meta_existing %>% group_by(biopsy_id) %>% tally() %>% filter(n>1) %>% select(biopsy_id)
  #singleton = rna_meta_existing %>% group_by(biopsy_id) %>% tally() %>% filter(n==1) %>% select(-n)

  #selected_duplicated = left_join(duplicated,rna_meta_existing,by="biopsy_id") %>%
  #  arrange(desc(protocol)) %>% arrange(biopsy_id) %>% group_by(biopsy_id) %>% filter(row_number()==1) %>% ungroup()
  #selected_singleton = left_join(singleton,rna_meta_existing,by="biopsy_id") %>% ungroup()
  #selected_all=rbind(selected_duplicated,selected_singleton)
  #put everything back together
  #at this point I have 1206 sample_ids
  ex_tidy = ex_tidy %>% dplyr::rename(mrna_sample_id = sample_id)

  ex_tidy = ex_tidy %>% dplyr::filter(mrna_sample_id %in% selected_libraries$sample_id)
  rna_meta = rna_meta %>% dplyr::select(sample_id,biopsy_id)

  ex_tidy_final = left_join(ex_tidy,rna_meta,by=c("mrna_sample_id"="sample_id"))
  #this still has the mrna sample ID. Need to add the tumour_sample_id from this biopsy (where available)
  genome_meta = get_gambl_metadata() %>%
    dplyr::select(biopsy_id,sample_id,patient_id,ffpe_or_frozen) %>%
    dplyr::rename("genome_sample_id" = "sample_id")
  #this gets the metadata in the same format but restricted to genome samples
  #REMOVE ANNOYING DUPLICATE GENOMES!
  duplicated_cases = genome_meta %>% group_by(biopsy_id) %>% tally() %>% dplyr::filter(n>1)

  selected_genomes = genome_meta %>%
    group_by(biopsy_id) %>%
    # Take the biopsy_id with the longest string length (e.g frozen)
    slice_max(str_length(ffpe_or_frozen), n = 1,with_ties=FALSE) %>%
    ungroup() %>% dplyr::select(-patient_id,-ffpe_or_frozen)

  #join to genome metadata based on biopsy_id (should be the same for RNA-seq and tumour genomes)
  ex_tidy_genome = left_join(ex_tidy_final,selected_genomes,by="biopsy_id")

  ex_tidy_genome = dplyr::select(ex_tidy_genome,Hugo_Symbol,mrna_sample_id,expression,biopsy_id,genome_sample_id) %>%
    dplyr::filter(!is.na(Hugo_Symbol))
  #write the data back out for use by others and loading into the database.
  tidy_expression_file = config::get("tidy_expression_file")

  #drop genes that are duplicated
  #table_out = pull(ex_tidy_genome,Hugo_Symbol) %>% table
  #badgenes = names(which(table_out>1803))
  #tidy_expression_data_nobadgenes = dplyr::filter(ex_tidy_genome,!Hugo_Symbol %in% badgenes)
  #full_expression_wider = dplyr::select(tidy_expression_data_nobadgenes,-mrna_sample_id,-genome_sample_id) %>%
  #    pivot_wider(names_from=Hugo_Symbol,values_from=expression)
  #transpose it
  #full_expression_wider = column_to_rownames(full_expression_wider,var="biopsy_id")
  #full_expression_wider_t = t(full_expression_wider)
  write_tsv(ex_tidy_genome,file=tidy_expression_file)


}


#' Title
#'
#' @param file_paths A vector of full file paths, e.g. the output of dir
#' @param tool_name The tool or pipeline that generated the files (should be the same for all)
#' @param output_type The file type to distinguish different output file types from the same pipeline (e.g. seg, maf, ploidy)
#' @param tool_version Optional: provide the version of the pipeline or tool
#' @param unix_group The unix group (should be the same for all)
#' @param sample_ids A vector of sample_id the same length and in the same order as the file paths
#' @param file_details_df Optionally supply the data frame directly instead (e.g. from find_files_extract_wildcards)
#'
#' @return Updates the database by appending to the gambl_files table. Use with caution!
#' @export
#'
#' @examples
assemble_file_details = function(file_details_df,file_paths,tool_name,unix_group,sample_ids,output_type="ploidy",is_production="yes"){
  #the gambl_files table contains
  #sample_id, unix_group, tool_name, tool_version, seq_type, genome_build, is_production, output_type, is_lifted_over, pairing_status
  # is_production, output_type, file_path, file_timestamp
  if(missing(file_details_df)){
  file_details_df = tibble(sample_id = sample_ids,unix_group=unix_group,tool_name=tool_name,output_type=output_type,is_production=is_production,
                    file_timestamp = lapply(file_paths,function(x){as.character(file.mtime(x))}),
                    file_path=file_paths) %>% unnest_longer(file_timestamp)
  }

  #date_info = file.mtime(file_path)
  con <- dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  dbWriteTable(con,"gambl_files",file_details_df,append=TRUE)
}


#' Use liftOver to convert a bedpe file between the two main genome builds (grch37/hg38)
#'
#' @param bedpe_file Either specify the path to a bedpe file
#' @param bedpe_df Or specify the bedpe data in a data frame
#' @param target_build Specify which build the data should be lifted to (must be one of hg19, grch37, hg38, grch38)
#'
#' @return Data frame containing original bedpe data with new coordinates
#' @export
#' @import tidyverse rtracklayer S4Vectors
#'
#' @examples
#' hg19_sv = get_manta_sv() %>% head(100)
#' hg38_sv = liftover_bedpe(bedpe_df=hg19_sv,target_build="hg38")
liftover_bedpe = function(bedpe_file,bedpe_df,target_build="grch37",col_names,col_types,standard_bed=FALSE){

  if(!missing(bedpe_file)){
    if(missing(col_names)){
      message("imposing column names")
      original_bedpe = read_tsv(bedpe_file,comment = "##",col_types="cddcddccccccccccccccccc")
    }else{
      message(paste("using column names",col_names,sep=": "))
      original_bedpe = read_tsv(bedpe_file,col_names=col_names,col_types=col_types)
    }
  }else{
    original_bedpe = bedpe_df

  }
  if(!standard_bed){
    colnames(original_bedpe)[1]="CHROM_A"
    original_bedpe = as.data.frame(original_bedpe)
    original_bedpe=original_bedpe %>% mutate_if(is.numeric, as.integer)
    #print(head(original_bedpe))
    if(!grepl("chr",original_bedpe$CHROM_A)){
      #add chr prefix
      original_bedpe = original_bedpe %>% mutate(CHROM_A = paste0("chr",CHROM_A)) %>%
      mutate(CHROM_B = paste0("chr",CHROM_B))
    }
    print(head(original_bedpe))
    char_vec = original_bedpe %>% tidyr::unite(united,sep="\t") %>% dplyr::pull(united)
    bedpe_obj <- rtracklayer::import(text=char_vec,format="bedpe")
    if(length(colnames(original_bedpe))>22){
      this_patient = colnames(original_bedpe)[23]
      this_normal = colnames(original_bedpe)[22]
    }else{
      this_patient = original_bedpe$tumour_sample_id
      this_normal = original_bedpe$normal_sample_id
    }
  }else{
    colnames(original_bedpe)[1]="chrom"
    if(!grepl("chr",original_bedpe$chrom)){
      original_bedpe = mutate(original_bedpe,chrom=paste0("chr",chrom))
    }
    char_vec = original_bedpe %>% tidyr::unite(united,sep="\t") %>% dplyr::pull(united)
    bedpe_obj <- rtracklayer::import(text=char_vec,format="bed")
  }
  if(target_build == "grch37" | target_build == "hg19"){
    chain = rtracklayer::import.chain(system.file("extdata","hg38ToHg19.over.chain",package="GAMBLR"))
  }else if(target_build == "grch38" | target_build == "hg38"){
    chain = rtracklayer::import.chain(system.file("extdata","hg19ToHg38.over.chain",package="GAMBLR"))
  }
  if(!standard_bed){
    colnames(original_bedpe)[1]="CHROM_A"
    original_columns = colnames(original_bedpe)

    first_sv_lifted = rtracklayer::liftOver(bedpe_obj@first,chain)
    second_sv_lifted = rtracklayer::liftOver(bedpe_obj@second,chain)
    no_problem = !((elementNROWS(first_sv_lifted) != 1) | (elementNROWS(second_sv_lifted) != 1))
    first_ok = subset(first_sv_lifted,no_problem)
    second_ok = subset(second_sv_lifted,no_problem)
    first_ok_df = rtracklayer::export(first_ok,format="bed") %>%
    read_tsv(col_names = c("CHROM_A","START_A","END_A","name_A","score_A","STRAND_A")) %>%
    dplyr::select(-score_A) %>% dplyr::select(-name_A)
    second_ok_df = rtracklayer::export(second_ok,format="bed") %>%
    read_tsv(col_names = c("CHROM_B","START_B","END_B","name_B","score_B","STRAND_B")) %>%
    dplyr::select(-score_B) %>% dplyr::select(-name_B)
    ok_bedpe = original_bedpe[no_problem,]
    kept_cols = ok_bedpe %>% dplyr::select(-c("CHROM_A","START_A","END_A","CHROM_B","START_B","END_B","STRAND_A","STRAND_B"))
    fully_lifted = bind_cols(first_ok_df,second_ok_df,kept_cols) %>% dplyr::select(all_of(original_columns))
    return(fully_lifted)
  }else{
    lifted = rtracklayer::liftOver(bedpe_obj,chain)
    no_problem = !((elementNROWS(lifted) != 1))
    first_ok = subset(lifted,no_problem)
    output = rtracklayer::export(first_ok,format="bed") %>%
      read_tsv(col_names = c("chrom","start","end","score","strand","nothing","s1","e1","junk","more","stuff","nada")) %>%
      dplyr::select("chrom","start","end")
    return(output)
  }
}





