#Global variable specifying what metadata columns are absolutely required
required_cols = c("sample_id","patient_id","pathology","seq_type","genome_build","pairing_status","Tumor_Sample_Barcode")

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




