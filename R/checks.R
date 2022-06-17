#Global variable specifying what metadata columns are absolutely required
required_cols = c("sample_id","patient_id","pathology","seq_type","genome_build","pairing_status","Tumor_Sample_Barcode")

#Helper functions not for export
grob_wildcards = function(wildcarded_string){
  wildcards = unlist(stringr::str_extract_all(wildcarded_string,"\\{[^\\{]+\\}"))
  wildcards = stringr::str_remove_all(wildcards,"\\{") %>%  stringr::str_remove_all(.,"\\}")
  return(wildcards)
}
get_template_wildcards = function(parent_key,template_key){
  if(missing(template_key)){
    wildcard_string = config::get(parent_key)
  }else{
    wildcard_string = config::get(paste0(parent_key,"_wildcards"))[template_key]
  }
  wildcards = stringr::str_split(wildcard_string,",")
  return(unlist(wildcards))
}

copy_no_clobber = function(from_file,to_file,force=FALSE){
  to_dir = dirname(to_file)
  suppressMessages(suppressWarnings(dir.create(to_dir,recursive = T)))
  print(paste("COPYING",from_file,"TO",to_file))
  if(force){
    file.copy(from_file,to_file)
  }
}

check_times = function(relative_paths,ssh_session,archive_mode=FALSE,force_backup=FALSE){
  
  local_base = base_path=config::get("project_base")
  remote_base = base_path=config::get("project_base",config="default")
  if(archive_mode){
    if(local_base == remote_base){
      message("checking against local archive")
    }else{
      message("Currently, this mode must be run on the GSC (not remotely)")
      return(NULL)
    }
    local_base = base_path=config::get("archive")
    
  }
  for(rel_f in relative_paths){
    local_f = paste0(local_base,rel_f)
    remote_f = paste0(remote_base,rel_f)
    print(rel_f)
    if(file.exists(local_f)){
      mtime = file.info(local_f)$mtime
      mtime = stringr::str_remove(mtime,"\\s\\d+:\\d+:\\d+")
      #print(mtime)
      
      
      #print(remote_f)
      if(!missing(ssh_session)){
        output = ssh::ssh_exec_internal(ssh_session,paste("stat -L ",remote_f,"| grep Modify"))$stdout
      
        output = rawToChar(output) %>% stringr::str_extract(.,"\\d+-\\d+-\\d+")
      }else{
        output = file.info(remote_f)$mtime %>% stringr::str_remove("\\s\\d+:\\d+:\\d+")
      }
      remote_time = lubridate::as_date(output)
      local_time = lubridate::as_date(mtime)
      agediff = lubridate::time_length(remote_time - local_time,unit="days")
      if(agediff>0){
        print(paste("Warning! Remote version is",agediff,"days newer than the local file. You probably need to update the following file:"))
        print(rel_f)
      }else{
        message("OK")
      }
    }else{
      if(archive_mode){
        message("local backup of this file doesn't seem to exist:")
        message(local_f)

        copy_no_clobber(remote_f,local_f,force_backup)

      }
    }
  }
}

check_file_details = function(relative_paths){
  not_found = c()
  base_path=config::get("project_base")
  for(relative_path in relative_paths){
    full_path = paste0(base_path,relative_path)
    message(paste("Looking for:",full_path))
    if(file.exists(full_path)){
      message("OK")
    }else{
      #message("Uh oh. This file cannot be found!")
      print(full_path)
      not_found=c(not_found,full_path)
      #print("-=10101010101=-")
    }
  }
  return(not_found)
}

check_gamblr_config = function(mode="default",compare_timestamps=FALSE,ssh_session,archive_mode=FALSE,force_backup=FALSE){
  files_to_check = c()
  #get all the wildcards we'll need
  seq_type = get_template_wildcards("seq_types")
  unix_group = get_template_wildcards("unix_groups")
  projection = get_template_wildcards("projections")
  
  #data frame for seq_type/projection expansion
  projection_expanded = tidyr::expand_grid(seq_type = get_template_wildcards("seq_types"),projection = get_template_wildcards("projections"))
  #resources section of config (only needs blacklist right now)
  blacklist_f = config::get("resources")$blacklist$template
  blacklist_f = mutate(projection_expanded,output=glue::glue(blacklist_f)) %>% pull(output)
  files_to_check = c(files_to_check,blacklist_f)
  
  merged_keys = names(config::get("results_merged"))
  #skip any file starting with "/"
  for (merge in merged_keys){
    merge_path = config::get("results_merged")[merge]
    if(!grepl("^/",merge_path)){
      files = unlist(merge_path)
      #print(names(files))
      
      for(f in files){
        #print(paste("CHECKING",f))
        if(grepl("\\{",f)){
          #print("contains wildcards, using all wildcards:")
          if(stringi::stri_count_fixed(f,"{")>1){
            print("Multiple wildcards!")
            wildcards = grob_wildcards(f)
            print(wildcards)
            if("projection" %in% wildcards & "seq_type" %in% wildcards){
              #use the same expansion approach we used above
              print(f)
              f = mutate(projection_expanded,output=glue::glue(f)) %>% pull(output)
              files_to_check = c(files_to_check,f)
            }else{
              warning("Unrecognized wildcards. Not sure how to handle this, using default approach")
              flavour =get_template_wildcards("results_merged",names(files))
              all_f = glue::glue(f)
              print(all_f)
              files_to_check = c(files_to_check,all_f)
            }
          }else{
            flavour =get_template_wildcards("results_merged",names(files))
            all_f = glue::glue(f)
            #print(all_f)
            files_to_check = c(files_to_check,all_f) 
          }
        }else{
          files_to_check = c(files_to_check,f)
        }
        
      }

    }
  }
  print(f)
  mia=check_file_details(files_to_check)
  l_missing = length(mia)
  if(l_missing){
    warning(paste("There were",l_missing,"files that cannot be found (see above). If this is unexpected, try to obtain them."))
    print("MISSING FILES:")
    print(mia)
  }
  if(compare_timestamps){
    check_times(files_to_check,ssh_session,archive_mode,force_backup)
  }
  print("DONE!")
}

#' Check GAMBL or other metadata for compatability with various features
#'
#' @param metadata_df Data frame output by get_gambl_metadata or some other source of metadata you plan to use
#' @param to_check Specify one of "uniqueness", "colours" or "completeness" or leave empty to check all
#' @param fix After identifying an issue, rerun this function with fix=TRUE to address errors (when possible). Currently this doesn't do anything. That's how I roll
#' @param show_details Set to TRUE if you want the gory details about issues that are identified
#'
#' @return
#' @export
#' @import dplyr
#'
#' @examples
check_gambl_metadata = function(metadata_df,to_check="all",show_details=FALSE,fix=FALSE){
  
  if(to_check == "all"){
    to_check = c("uniqueness","colours","completeness")
  }
  if("completeness" %in% to_check){
    if(any(!required_cols %in% colnames(metadata_df))){
      missing = required_cols[!required_cols %in% colnames(metadata_df)]
      message("FAIL! MISSING metadata columns:")
      print(paste(missing,sep=" "))
      return()
    }
  }
  if("uniqueness" %in% to_check){
    # check that there are no duplicate sample_id and warn the user of any violations thereof
    if(any(duplicated(metadata_df$sample_id))){
      numdup = sum(duplicated(metadata_df$sample_id))
      warning(paste("Some",numdup,"values in your sample_id column are duplicates"))
    }
    if(any(duplicated(metadata_df$Tumor_Sample_Barcode))){
      numdup = sum(duplicated(metadata_df$Tumor_Sample_Barcode))
      warning(paste("Some",numdup,"values in your sample_id column are duplicates"))
    }
    message("PASSED uniqueness test for sample_id")
    
  }
  if("colours" %in% to_check){
    # confirm that we have colours for all values in all columns that map to a colour when we run map_metadata_to_colours
    # For this to work, they'll need an alias in the global variable colour_aliases. Users can add their own aliases if they want them handled
    total=nrow(metadata_df)
    alias_names = names(colour_aliases)
    alias_in_meta = alias_names[alias_names %in% colnames(metadata_df)]
    #print(paste("will check for colour mapping of values in",alias_in_meta))
    for(alias in alias_in_meta){
      mapped = data.frame(map_metadata_to_colours(alias,metadata_df,as_vector=F)) %>% 
        rename("colour"=alias) %>% 
        rownames_to_column(var=alias)
      
      this_col = pull(metadata_df,alias)
      
      if(class(this_col)=="factor"){
        metadata_df = mutate(metadata_df,{{alias}} := as.character(.data[[alias]]))
      }
      tabular = group_by(metadata_df,!!!syms(alias)) %>% 
        tally() %>%
        mutate({{alias}} := replace_na(.data[[alias]],"NA"))
      mapped = right_join(mapped,tabular,by=alias) %>% dplyr::filter(!is.na(n))
      
      sum_na = dplyr::filter(mapped,is.na(colour)) %>% pull(n) %>% sum()
      percent=round(100 * sum_na/total,3)
      if(sum_na ==0){
        message(paste(alias,"OK"))
      }else{
        message("possible problem!")
        print(paste(percent,"% of values not assigned to an available colour."))
        if(show_details){
          bad = dplyr::filter(mapped,is.na(colour))
          print(bad)
        }
      }
    }
  }
  if(fix){
    message("Nothing fixed. I'm waiting for someone to implement some code to fix common issues.")
    message("maybe try something like this? GAMBLR::tidy_lymphgen(metadata,  lymphgen_with_cnv, lymphgen_with_cnv_tidy, relevel = TRUE)")
  }
}

get_runs_table = function(seq_type_filter="genome"){
  t_meta = get_gambl_metadata(tissue_status_filter = c("tumour"),seq_type_filter=seq_type_filter) %>% 
    dplyr::select(sample_id,patient_id,seq_type,genome_build,pairing_status) %>% rename("tumour_sample_id"="sample_id")
  n_meta = get_gambl_metadata(tissue_status_filter = c("normal"),seq_type_filter=seq_type_filter) %>%
    dplyr::select(sample_id,patient_id,seq_type,genome_build) %>% rename("normal_sample_id"="sample_id")
}

check_expected_outputs = function(tool_name="battenberg",seq_type_filter="genome"){
  all_meta = get_gambl_metadata(seq_type_filter=seq_type_filter)
  #drop irrelevant rows of the metadata based on the scope of the tool etc
  if(tool_name=="battenberg"){
    template_path = config::get("results_flatfiles")$cnv$battenberg
    extra_wildcards = config::get("results_flatfiles")$cnv$battenberg_wildcards
    #in the current setup, this drops unmatched samples (could be hard-coded but using the config is more transparent)
    relevant_metadata = dplyr::filter(all_meta,base::get(names(extra_wildcards)[1]) == unname(extra_wildcards[1])) 
    relevant_metadata = dplyr::filter(relevant_metadata,seq_type == seq_type_filter)
  }
  w = grob_wildcards(template_path)

  projection = get_template_wildcards("projections")
  #use the unix group and seq_type from the actual metadata and all available projections
  seq_type = seq_type_filter
  template_path = glue::glue()
}


