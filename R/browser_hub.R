
deblacklist_regions = function(regions_bed,projection="grch37"){
  if(projection=="grch37"){
    blacklist_file = "/Users/rmorin/git/LLMPP/resources/reference/encode/hg19-blacklist.v2.bed"
  }
  blacklist_df = read_tsv(blacklist_file,col_names = c("chr","start","end","name"))
  blacklist_dt = as.data.table(blacklist_df)
  colnames(regions_bed)[c(1:3)]=c("chr","start","end")
  if(any(!grepl("chr",regions_bed[,1]))){
    message("adding chr prefix")
    regions_bed = mutate(regions_bed,chr = as.character(chr))
    regions_bed = mutate(regions_bed,chr = paste0("chr",chr))
  }
  regions_dt = as.data.table(regions_bed)

 
  #annotate bins that overlap blacklisted regions
  setkey(blacklist_dt, chr,start,end)
  setkey(regions_dt, chr,start,end)

  annotated = foverlaps(blacklist_dt, regions_dt,by.x =c("chr","start","end"),by.y=c("chr","start","end"),mult = "all",type = "any") %>% dplyr::filter(!is.na(start))
  to_keep = dplyr::filter(regions_bed,!start %in% annotated$start)
  return(to_keep)
}

load_combine_cached_bins = function(function_params=list(bin_size=15000),
  additional_details=list(foreground="DLBCL_FL",background="CLL_MM")){

      cache_file_name = paste0(config::get("repo_base"),"cached_results/get_ssm_by_region")
      for (param in names(function_params)[order(names(function_params))]){
        cache_file_name = paste0(cache_file_name,"--",param,"-",function_params[[param]])
      }
      cache_file_name = paste0(cache_file_name,"*")
      for (detail in names(additional_details)){
        cache_file_name = paste0(cache_file_name,"--",detail,"-",additional_details[[detail]])
      }
      cache_file_name = paste0(cache_file_name,".tsv")
     
    cached = Sys.glob(paths = cache_file_name )
  if(!length(cached)==48){
    stop(paste("error: cannot find 48 files using this pattern",cache_file_name))
  }
  binned_list = list()
  for(to_add in cached){
    this_df = read_tsv(to_add,col_types="ciicnnnn")
    binned_list[[to_add]]=this_df
  }
  return(bind_rows(binned_list))
}

bins_to_bedgraph = function(bin_df,min_value = 3,filename = "test.bed"){
  bed_cols = dplyr::select(bin_df,1,2,3,smoothed_ratio)
  colnames(bed_cols) = c("chr","start","end","value")
  if(any(!grepl("chr",bed_cols[,1]))){
    message("adding chr prefix")
    bed_cols = mutate(bed_cols,chr = paste0("chr",chr))
  }
  bed_cols = dplyr::filter(bed_cols,value > min_value)
  this = dplyr::filter(chr1p,bin_start==27985000) 
  bed_cols = mutate(bed_cols,end=format(end, scientific=F),start=format(start, scientific=F))
  write.table(bed_cols,row.names=F,col.names=F,quote=F,sep="\t",file=filename)
  
}

build_browser_lymphomagene = function(browser_base_path="/Users/rmorin/git/LLMPP/",
                              projection="grch37",
                              hub_path="hubs/ashm/",
                              ucsc_tools_path="/Users/rmorin/miniconda3/envs/ucsc/bin/"){
  if(projection=="grch37"){
    genome_build = "hg19"
    # step 0: start writing to the hub.txt file (clobbering existing file)
    hubfile = paste0(browser_base_path,hub_path,"hub.txt")
    cat("hub aSHM sites in hg19\nshortLabel aSHM_hg19\nlongLabel curated aSHM sites in hg19\nuseOneFile on\nemail rdmorin@sfu.ca\n\ngenome hg19\n\n",file=hubfile)
    # step 1: format the curated SHM regions file
    ashm_file = paste0(browser_base_path,"resources/curated/somatic_hypermutation_locations_GRCh37.txt")
    ashm_bed = gsub(".txt",".bed",ashm_file)
    #strip additional details from file and ensure it's sorted
    ashm_details = read_tsv(ashm_file) %>% select(1,2,3) %>% 
      arrange(1,2)
    message(paste("overwriting:",ashm_bed))
    write.table(ashm_details,file=ashm_bed,sep="\t",col.names  = F,row.names=F,quote=F)
    #convert to a bigBed file in browser directory
    bigbed = paste0(browser_base_path,hub_path,"somatic_hypermutation_locations_GRCh37.bb")
    ref_sizes = paste0(browser_base_path,"resources/reference/ucsc/",genome_build,".chrom.sizes")
    conversion = paste0(ucsc_tools_path,"bedToBigBed"," ",ashm_bed," ",ref_sizes," ",bigbed)
    system(conversion)
    cat("track aSHM_regions\nshortLabel aSHM regions\nlongLabel curated aSHM sites\nvisibility full\ntype bigBed\nbigDataUrl https://github.com/morinlab/LLMPP/blob/main/hubs/ashm/somatic_hypermutation_locations_GRCh37.bb?raw=true\n\n",file=hubfile,append=TRUE)
    
    # step 2: obtain all mutations in these regions 
    ashm_regions_bed = read_tsv(ashm_file)
    genome_ssm = get_ssm_by_regions(regions_bed = ashm_regions_bed,streamlined = F,seq_type = "genome",projection = projection)
    ashm_full_bed = paste0(browser_base_path,hub_path,"ashm_",genome_build,".bb")
    maf_to_custom_track(genome_ssm,
                        these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome"),
                        as_bigbed = T,output_file = ashm_full_bed
    )
    # step 2.1: update hub file with information for new track
    cat("track mutations_aSHM\n",file=hubfile,append=TRUE)
    cat("shortLabel aSHM mutations\n",file=hubfile,append=TRUE)
    cat("longLabel mutations at aSHM sites\n",file=hubfile,append=TRUE)
    cat("visibility squish\n",file=hubfile,append=TRUE)
    cat("type bigBed 9\nitemRgb on\n",file=hubfile,append=TRUE)
    cat("bigDataUrl https://github.com/morinlab/LLMPP/blob/main/hubs/ashm/ashm_hg19.bb?raw=true\n",append=TRUE,file=hubfile)
    
    # step 3: annotate recurrent mutations for lollipop
    
    coding_genome_variants = get_coding_ssm(seq_type = "genome")
    
  }  
}

build_browser_ashm = function(browser_base_path="/Users/rmorin/git/LLMPP/",
                              projection="grch37",
                              hub_path="hubs/ashm/",
                              ucsc_tools_path="/Users/rmorin/miniconda3/envs/ucsc/bin/"){
  if(projection=="grch37"){
    genome_build = "hg19"
    # step 0: start writing to the hub.txt file (clobbering existing file)
    hubfile = paste0(browser_base_path,hub_path,"hub.txt")
    cat("hub aSHM sites in hg19\nshortLabel aSHM_hg19\nlongLabel curated aSHM sites in hg19\nuseOneFile on\nemail rdmorin@sfu.ca\n\ngenome hg19\n\n",file=hubfile)
    # step 1: format the curated SHM regions file
    ashm_file = paste0(browser_base_path,"resources/curated/somatic_hypermutation_locations_GRCh37.txt")
    ashm_bed = gsub(".txt",".bed",ashm_file)
    #strip additional details from file and ensure it's sorted
    ashm_details = read_tsv(ashm_file) %>% select(1,2,3) %>% 
      arrange(1,2)
    message(paste("overwriting:",ashm_bed))
    write.table(ashm_details,file=ashm_bed,sep="\t",col.names  = F,row.names=F,quote=F)
    #convert to a bigBed file in browser directory
    bigbed = paste0(browser_base_path,hub_path,"somatic_hypermutation_locations_GRCh37.bb")
    ref_sizes = paste0(browser_base_path,"resources/reference/ucsc/",genome_build,".chrom.sizes")
    conversion = paste0(ucsc_tools_path,"bedToBigBed"," ",ashm_bed," ",ref_sizes," ",bigbed)
    system(conversion)
    cat("track aSHM_regions\nshortLabel aSHM regions\nlongLabel curated aSHM sites\nvisibility full\ntype bigBed\nbigDataUrl https://github.com/morinlab/LLMPP/blob/main/hubs/ashm/somatic_hypermutation_locations_GRCh37.bb?raw=true\n\n",file=hubfile,append=TRUE)

    # step 2: obtain all mutations in these regions 
    ashm_regions_bed = read_tsv(ashm_file)
    genome_ssm = get_ssm_by_regions(regions_bed = ashm_regions_bed,streamlined = F,seq_type = "genome",projection = projection)
    ashm_full_bed = paste0(browser_base_path,hub_path,"ashm_",genome_build,".bb")
    maf_to_custom_track(genome_ssm,
                        these_samples_metadata = get_gambl_metadata(seq_type_filter = "genome"),
                        as_bigbed = T,output_file = ashm_full_bed
                        )
    # step 2.1: update hub file with information for new track
    cat("track mutations_aSHM\n",file=hubfile,append=TRUE)
    cat("shortLabel aSHM mutations\n",file=hubfile,append=TRUE)
    cat("longLabel mutations at aSHM sites\n",file=hubfile,append=TRUE)
    cat("visibility squish\n",file=hubfile,append=TRUE)
    cat("type bigBed 9\nitemRgb on\n",file=hubfile,append=TRUE)
    cat("bigDataUrl https://github.com/morinlab/LLMPP/blob/main/hubs/ashm/ashm_hg19.bb?raw=true\n",append=TRUE,file=hubfile)
    
    # step 3: annotate recurrent mutations for lollipop
    
    coding_genome_variants = get_coding_ssm(seq_type = "genome")
    
  }  
}
