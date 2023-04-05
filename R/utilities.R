#global environment
coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")
rainfall_conv = c("T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G", "InDel")
names(rainfall_conv) = c('A>G', 'T>C', 'C>T', 'G>A', 'A>T', 'T>A', 'A>C', 'T>G', 'C>A', 'G>T', 'C>G', 'G>C', 'InDel')
ssh_session <<- NULL


compare_coding_mutation_pattern = function(maf_df1,maf_df2,gene){
  if(missing(maf_df1) | missing(maf_df2)){
    stop("must provide two data frames containing mutations you would like to compare")
  }
  if(missing(gene)){
    stop("Must provide the Hugo_Symbol of a single gene that is present in both maf files")
  }
  missense_positions1 = dplyr::filter(maf_df1,Hugo_Symbol==gene,!Variant_Classification %in% c("Silent","Splice_Site","Splice_Region"),Variant_Type=="SNP") %>%
    pull(HGVSp_Short) %>% str_remove_all("p.\\w") %>% str_extract("\\d+") %>% as.numeric()
  missense_positions2 = dplyr::filter(maf_df2,Hugo_Symbol==gene,!Variant_Classification %in% c("Silent","Splice_Site","Splice_Region"),Variant_Type=="SNP") %>%
    pull(HGVSp_Short) %>% str_remove_all("p.\\w") %>% str_extract("\\d+") %>% as.numeric()
 if(length(missense_positions1)==0 | length(missense_positions2)==0 ){
   message(paste("no mutations for",gene,"in one or both data sets"))
   return(list(kl=15))
 }
  #generate range of amino acids based on what we can infer from the MAF (not ideal)
  max_pos = max(c(missense_positions1,missense_positions2))
  full_df = data.frame(position=c(1:max_pos))
  df1 = data.frame(position=missense_positions1) %>% group_by(position) %>% tally() %>% rename("group1"="n")
  df2 = data.frame(position=missense_positions2) %>% group_by(position) %>% tally() %>% rename("group2"="n")
  full_df = left_join(full_df,df1) %>% mutate(group1=ifelse(is.na(group1),0,group1))
  full_df = left_join(full_df,df2) %>% mutate(group2=ifelse(is.na(group2),0,group2))
  # convert to the format needed by KL
  all_counts = dplyr::select(full_df,-position) %>% t()
  all_counts[1,]=all_counts[1,]/sum(all_counts[1,])
  all_counts[2,]=all_counts[2,]/sum(all_counts[2,])
  kl_out = KL(all_counts)
  return(list(df=full_df,kl=unname(kl_out)))
}

#' @title Write Sample Set Hash
#'
#' @description Update or create a file to track unique identifiers for sample sets in GAMBL
#'
#' @details Run this function with `update = TRUE` (default) to use an existing sample table.
#' If this table does not exist, perhaps you need to pull from the master branch.
#' If this function is run with the default for `update`, the user must also provide the new sample sets with the `new_sample_sets_df`.
#'
#' @param update Leave as TRUE for default functionality (i.e. updating the existing table). If the table doesn't exist you probably need to pull from Master.
#' @param new_sample_sets_df Data frame of all existing and new sample sets. Required when running in default update mode.
#'
#' @import dplyr readr
#' @export
#'
write_sample_set_hash = function(update = TRUE,
                                 new_sample_sets_df){

  sample_sets_file = paste0(check_config_value(config::get("repo_base")), check_config_value(config::get("sample_sets")$default))
  md5_file = paste0(check_config_value(config::get("repo_base")), check_config_value(config::get("sample_sets")$hashes))

  if(update){
    # load the existing file and update it using the contents of sample_sets_df as well as checking for consistency for existing sample sets
    if(missing(new_sample_sets_df)){
      stop("You must provide a data frame containing all the sample sets to update the digests")
    }
    original_digests = suppressMessages(read_tsv(md5_file))

    set_names = dplyr::select(new_sample_sets_df,-sample_id) %>% colnames()
    #only compare for sample sets that we have in the current file
    md5_values= c()
    for(set_name in set_names){

      this_md5 = get_samples_md5_hash(sample_set_name=set_name,sample_sets_df = new_sample_sets_df)
      md5_values=c(md5_values,this_md5)

    }
    all_md5 = data.frame(sample_set = set_names,new_md5_digest=md5_values)
    oldnew = right_join(original_digests,all_md5,by="sample_set")
    #check the rows where md5_digest is not NA (i.e. rows that were there before)
    to_check = dplyr::filter(oldnew,!is.na(md5_digest))
    if(any(to_check$md5_digest != to_check$new_md5_digest)){
      problems = dplyr::filter(to_check,md5_digest!=new_md5_digest)
      print(problems)
      stop("some md5 digests do not match. Have these sample sets changed???")

    }

  }else{
    # just create a file that records the md5 digests for existing sample sets
    sample_sets = suppressMessages(read_tsv(sample_sets_file))
    set_names = select(sample_sets,-sample_id) %>% colnames()
    message(paste("Will log digests for",length(set_names),"sample sets"))
    md5_values= c()
    for(set_name in set_names){
      this_md5 = get_samples_md5_hash(sample_set_name=set_name)
      md5_values=c(md5_values,this_md5)
    }
    all_md5 = data.frame(sample_set = set_names,md5_digest=md5_values)
    write_tsv(all_md5,file=md5_file)
  }
}


#' @title Generate md5 Hash For Samples
#'
#' @description Generate an md5 hash for a set of samples to help ensure reproducibility
#'
#' @details This function can accept a wide range of formatted sample IDs to create an md5 hash.
#' For example, if the user is working with an already subset metadata table (with sample IDs of interest),
#' The user can give this table to the function with `these_sampels_metadata`.
#' As an alternative, sample IDs can also be provided as a vector of characters with `these_samples` parameter.
#' Another option is to use defined sample sets (GAMBL) with `sample_set_name`.
#' As a final option, the user can also provide a data frame with samples IDs instead of loading them from the GAMBL repo,
#' This is achieved with calling the `sample_sets_df` parameter.
#'
#'
#' @param these_samples_metadata Optionally provide a metadata table or any data frame with a column named sample_id that has been subset to the samples you're working with.
#' @param these_samples Optionally provide a vector of sample_id you are working with.
#' @param sample_set_name Optionally provide the name of a sample set in GAMBL and the function will load the samples from that set and provide the hash.
#' @param sample_sets_df Optionally provide a data frame of the sample sets instead of relying on/loading the local file from the GAMBL repo.
#'
#' @return The md5 hash of the ordered set of sample_id.
#'
#' @import digest dplyr readr
#' @export
#'
get_samples_md5_hash = function(these_samples_metadata,
                                these_samples,
                                sample_set_name,
                                sample_sets_df){

  if(!missing(these_samples_metadata)){
    collapsed = dplyr::select(these_samples_metadata,sample_id) %>%
      arrange() %>%
      pull() %>% paste(.,collapse=",")

      digested = digest::digest(collapsed,serialize = FALSE)
  }else if(!missing(these_samples)){
    digested = digest::digest(paste(these_samples[order(these_samples)],collapse=","),serialize=FALSE)
  }else if(!missing(sample_set_name)){
    #load the sample set table and pull the samples based on its contents and the name provided
    sample_sets_file = paste0(check_config_value(config::get("repo_base")), check_config_value(config::get("sample_sets")$default))
    if(missing(sample_sets_df)){
      sample_sets = suppressMessages(read_tsv(sample_sets_file))
    }else{
      sample_sets = sample_sets_df
    }
    setname = as.symbol(sample_set_name)
    collapsed = dplyr::filter(sample_sets,!!setname==1) %>%
      dplyr::select(sample_id) %>%
      arrange() %>% pull() %>% paste(.,collapse=",")
    digested = digest::digest(collapsed,serialize=FALSE)
  }
  return(digested)
}


cache_output = function(result_df,
                        function_name,
                        clobber_mode = F,
                        get_existing = F,
                        function_params = list(region = "chr3:98300000-198022430", bin_size=2000, seq_type="genome"),
                        additional_details = list(foreground = "DLBCL_FL_BL", background = "CLL_MM_MCL")){

  cache_file_name = paste0(check_config_value(config::get("repo_base")),"cached_results/", function_name)
  for (param in names(function_params)[order(names(function_params))]){
    cache_file_name = paste0(cache_file_name,"--",param,"-",function_params[[param]])
  }
  for (detail in names(additional_details)){
    cache_file_name = paste0(cache_file_name,"--",detail,"-",additional_details[[detail]])
  }
  cache_file_name = paste0(cache_file_name,".tsv")
  if(file.exists(cache_file_name)){
    if(get_existing){
      result_df = suppressMessages(read_tsv(cache_file_name))
      return(result_df)
    }
    if(!clobber_mode){
      warning(paste("file",cache_file_name,"exists!"))
      stop("Will not overwrite unless you rerun this in clobber_mode = TRUE")
    }
  }else{
    if(get_existing){
      stop(paste("cannot find cached result for this parameter combination",cache_file_name))
    }
  }

  message(paste("creating/overwriting",cache_file_name))
  write_tsv(result_df,file=cache_file_name)
}

#' @title Count SSM In A Region
#'
#' @description Count the variants in a region with a variety of filtering options.
#'
#' @details This function internally calls [GAMBLR::get_ssm_by_region] thus, the parameters available to this function are arguments that are being passed to the internal call.
#' For more details on how these parameters can be used, refer to [GAMBLR::get_ssm_by_region].
#'
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param start Query start coordinate of the range you are restricting to.
#' @param end Query end coordinate of the range you are restricting to.
#' @param these_samples_metadata A metadata table subset to the sample IDs of interest. If not provided, the function will call `get_gambl_metadata` and regions will be returned for all samples in the metadata.
#' @param all_mutations_in_these_regions If you are calling this function many times (e.g. bins spanning a larger region), to save a ton of time you are strongly encouraged to provide the output of `get_ssm_by_region` on the entire region of interest and passing it to this function
#' @param count_by Defaults to counting all variants. Specify 'sample_id' if you want to collapse and count only one per sample
#' @param seq_type The seq_type you want back, default is genome.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' library(dplyr)
#' 
#' #define a region.
#' my_region = gene_to_region(gene_symbol = "MYC", 
#'                            return_as = "region")
#'
#' #subset metadata.
#' my_metadata = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL")
#'
#' #count SSMs for the selected sample subset and defined region.
#' fl_ssm_counts_myc = count_ssm_by_region(region = my_region,
#'                                        these_samples_metadata = my_metadata)
#'
count_ssm_by_region = function(region,
                               chromosome,
                               start,
                               end,
                               all_mutations_in_these_regions,
                               these_samples_metadata,
                               count_by,
                               seq_type = "genome"){

  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter = seq_type)
  }
  if(!missing(all_mutations_in_these_regions)){
    # function was provided the mutations already so we just need to subset it to the region of interest

    region_muts = dplyr::filter(all_mutations_in_these_regions,Start_Position >= start, Start_Position < end)
  }else if(missing(region)){
    region_muts = get_ssm_by_region(chromosome=chromosome,qstart=start,qend=end,streamlined = TRUE)
  }else{
    region_muts = get_ssm_by_region(region=region,streamlined = TRUE)
  }
  keep_muts = dplyr::filter(region_muts,Tumor_Sample_Barcode %in% these_samples_metadata$Tumor_Sample_Barcode)
  if(missing(count_by)){
    #count everything even if some mutations are from the same patient
    return(nrow(keep_muts))
  }else if(count_by == "sample_id"){
    return(length(unique(keep_muts$Tumor_Sample_Barcode)))
  }else{
    print("Not sure what to count")
  }
}

#' @title Region To Bins.
#'
#' @description Split a contiguous genomic region on a chromosome into non-overlapping bins
#'
#' @details This function takes genomic coordinates with the `chromosome`, `start`, and `end` parameters.
#' Lastly, the user can also specify the bin size with `bin_size`. Default is 20000
#'
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param start Query start coordinate of the range you are restricting to.
#' @param end Query end coordinate of the range you are restricting to.
#' @param bin_size The size of the bins, default is 2000.
#'
#' @return Data frame describing the bins various ways.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' chr8q_bins = region_to_bins(chromosome = "8",
#'                             start = 48100000,
#'                             end = 146364022,
#'                             bin_size = 20000)
#' 
region_to_bins = function(chromosome,
                 start,
                 end,
                 bin_size = 2000){

  if(missing(chromosome)){
    stop("Please provide a chromosome...")
  }

  if(missing(start)){
    stop("Please provide the start coordinates...")
  }

  if(missing(end)){
    stop("Please provide the end coordinates...")
  }

  bin_df = data.frame(bin_chr = chromosome, bin_start = seq(start, end, bin_size))
  bin_df = mutate(bin_df,bin_end = bin_start+ bin_size) %>%
    dplyr::filter(bin_end<=end) %>%
    mutate(region = paste0(bin_chr,":",bin_start,"-",bin_end))

  return(bin_df)
}


#' @title Get SSH Session.
#'
#' @description Create an ssh session to the GSC (requires active VPN connection)
#'
#' @details Using the ssh R package to create an ssh session.
#'
#' @param host Default is "gphost01.bcgsc.ca".
#'
#' @return An external pointer of class 'ssh_session'
#'
#' @import ssh
#' @export
#'
#' @examples
#' my_session = get_ssh_session()
#'
get_ssh_session = function(host="gphost01.bcgsc.ca"){

  if(!is.null(config::get("host"))){
    host = config::get("host")
  }

  if (!requireNamespace("ssh", quietly = TRUE)) {
    warning("The ssh package must be installed to use this functionality")
    #Either exit or do something that does not require ssh
    return(NULL)
  }
  message("you should also run this command to ensure the ssh library is loaded:\nlibrary(ssh)")
  session = ssh::ssh_connect(host=host)
  return(session)
}


#' @title Gene To Region.
#'
#' @description Return coordinates for a given gene or a set of genes.
#'
#' @details This function takes one or multiple gene names, either as hugo symbols or Ensembl IDs
#' and returns the coordinates for the selected genes in multiple formats (`return_as`).
#' The possible return formats are; bed, data frame and in "region" format (chr:start-end).
#' For returning genes residing in specific regions, see [GAMBLR::region_to_gene].
#'
#' @param gene_symbol A vector of one or more gene symbols.
#' @param ensembl_id A vector of one or more Ensembl IDs.
#' @param genome_build Reference genome build, default is grch37.
#' @param return_as Specify the type of return. Default is region (chr:start-end), other acceptable arguments are "bed" and "df".
#'
#' @return A data frame, or a string with region(s) for the provided gene(s).
#'
#' @import dplyr
#' @export
#'
#' @examples
#' bcl2_region = gene_to_region(gene_symbol = "BCL2",
#'                              genome_build = "grch37")
#' 
#' bcl2_region = gene_to_region(ensembl_id = "ENSG00000171791",
#'                              genome_build = "grch37")
#'
gene_to_region = function(gene_symbol,
                          ensembl_id,
                          genome_build = "grch37",
                          return_as = "region"){

  #set mart based on selected genome projection
  if(genome_build == "grch37"){
    gene_coordinates = grch37_gene_coordinates
    chr_select = paste0(c(c(1:22),"X","Y"))
  }else if(genome_build == "hg38"){
    gene_coordinates = hg38_gene_coordinates
    chr_select = paste0("chr", c(c(1:22),"X","Y"))
  }

  #filter on gene_symbol/ensembl_id
  if(!missing(gene_symbol) && missing(ensembl_id)){
    gene_coordinates = dplyr::filter(gene_coordinates, hugo_symbol %in% gene_symbol)
  }

  if(missing(gene_symbol) && !missing(ensembl_id)){
    gene_coordinates = dplyr::filter(gene_coordinates, ensembl_gene_id %in% ensembl_id)
  }

  region = dplyr::select(gene_coordinates, chromosome, start, end, gene_name, hugo_symbol, ensembl_gene_id) %>%
    as.data.frame() %>%
    dplyr::arrange(chromosome, start) %>%
    dplyr::filter(chromosome %in% chr_select) %>%
    mutate_all(na_if,"") %>%
    distinct(.keep_all = TRUE)

  if(return_as == "bed"){
    #return one-row data frame with first 4 standard BED columns. TODO: Ideally also include strand if we have access to it in the initial data frame
    region = dplyr::select(region, chromosome, start, end, hugo_symbol)

  }else if(return_as == "df"){
    region = region

  }else{
    #default: return in chr:start-end format
    region = paste0(region$chromosome, ":", region$start, "-", region$end)
  }

  if(return_as %in% c("bed", "df")){
    if(!missing(gene_symbol)){
      message(paste0(nrow(region), " region(s) returned for ", length(gene_symbol), " gene(s)"))
    }

    if(!missing(ensembl_id)){
      message(paste0(nrow(region), " region(s) returned for ", length(ensembl_id), " gene(s)"))
    }
  }else{
    if(!missing(gene_symbol)){
      message(paste0(length(region), " region(s) returned for ", length(gene_symbol), " gene(s)"))
    }

    if(!missing(ensembl_id)){
      message(paste0(length(region), " region(s) returned for ", length(ensembl_id), " gene(s)"))
    }
  }
  return(region)
}


#' @title Region To Gene.
#'
#' @description Return genes residing in defined region(s).
#'
#' @details This function takes a region as a vector of characters, or a data frame with of regions (e.g output from `gene_to_region(return_as="df")`).
#' and returns the genes residing withing the specified region. For the other way around (i.e gene to regions, refer to [GAMBLR::gene_to_region]).
#'
#' @param region Regions to intersect genes with, this can be either a data frame with regions sorted in the following columns; chromosome, start, end. Or it can be a character vector in "region" format, i.e chromosome:start-end.
#' @param gene_format Parameter for specifying the format of returned genes, default is "hugo", other acceptable inputs are "ensembl".
#' @param genome_build Reference genome build.
#'
#' @return A data frame with gene(s) that are residing in the specified region(s).
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr
#' @export
#'
#' @examples
#' myc_region = gene_to_region(gene_symbol = "MYC",
#'                             genome_build = "grch37",
#'                             return_as = "df")
#' 
#' region = region_to_gene(region = myc_region,
#'                         gene_format = "hugo",
#'                         genome_build = "grch37")
#'
region_to_gene = function(region,
                          gene_format = "hugo",
                          genome_build = "grch37"){

  #set mart based on selected genome projection
  if(genome_build == "grch37"){
    gene_list = grch37_gene_coordinates
  }else if(genome_build == "hg38"){
    gene_list = hg38_gene_coordinates
  }

  #rename columns to match downstream formats
  colnames(gene_list)[1] = "ensembl_gene_id"
  colnames(gene_list)[2] = "chromosome"
  colnames(gene_list)[3] = "start"
  colnames(gene_list)[4] = "end"
  colnames(gene_list)[5] = "gene_name"
  colnames(gene_list)[6] = "hugo_symbol"

  gene_list = as.data.frame(gene_list)

  if(is.data.frame(region)){
    region_table = as.data.table(region)
  }else if(is.character(region)){
    split_chunks = unlist(strsplit(region, ":"))
    split_chunks = unlist(strsplit(split_chunks, "-"))
    chromosome = split_chunks[1]
    start = split_chunks[2]
    end = split_chunks[3]
    region = cbind(chromosome, start, end) %>%
      as.data.frame()

    region_table = as.data.table(region)

    region_table$chromosome = as.character(region_table$chromosome)
    region_table$start = as.double(region_table$start)
    region_table$end = as.double(region_table$end)
  }

  #transform regions to data tables
  gene_table = as.data.table(gene_list)

  #set keys
  data.table::setkey(region_table, chromosome, start, end)
  data.table::setkey(gene_table, chromosome, start, end)

  #intersect regions
  intersect = data.table::foverlaps(region_table, gene_table, nomatch = 0)
  colnames(intersect)[7] = "region_start"
  colnames(intersect)[8] = "region_end"

  #transform object to data frame
  inter_df = as.data.frame(intersect)

  #organize columns to match the expected format
  if(gene_format == "hugo"){
    genes = dplyr::select(inter_df, chromosome, start, end, hugo_symbol, region_start, region_end)
  }else if(gene_format == "ensembl"){
    genes = dplyr::select(inter_df, chromosome, start, end, ensembl_gene_id, region_start, region_end)
  }

  #paste chr in chromosome column, if not there
  if(!str_detect(genes$chromosome[1], "chr")){
    genes = mutate(genes, chromosome = paste0("chr", chromosome))}

  genes = as.data.frame(genes) %>%
    dplyr::arrange(chromosome, start) %>%
    distinct(.keep_all = TRUE)

  message(paste0(nrow(genes), " gene(s) returned for ", nrow(region), " region(s)"))

  return(genes)
}


#' @title Compare Mutation Flavour.
#'
#' @description Get a MAF that is just the variants unique to one of two flavours of variant calls available.
#'
#' @details Subset a MAF to only have variants that are unique to one flavour (specified with `flavour1`).
#' This function is currently not exported, since there is only one flavour available at the moment (see docs for [GAMBLR::get_ssm_by_sample]).
#'
#' @param these_sample_ids A vector of sample IDs to be included.
#' @param flavour1 First flavour of variant calls, to be returned as unique if not present in flavour2. Default is "clustered".
#' @param flavour2 Second flavour of variant calls.
#'
#' @return a list with MAFs that are only present in flavour1.
#'
#' @noRd
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#'
compare_mutation_flavour = function(these_sample_ids,
                                    flavour1 = "clustered",
                                    flavour2 = ""){

  these_dfs = list()
  for(this_sample_id in these_sample_ids){
    message(this_sample_id)
    maf1 = get_ssm_by_sample(this_sample_id, flavour = flavour1, this_seq_type = seq_type)
    maf2  = get_ssm_by_sample(this_sample_id, flavour = flavour2)
    maf1_only = intersect_maf(maf1, maf2)
    these_dfs[[this_sample_id]] = maf1_only
  }
  this_maf = rbindlist(these_dfs, use.names = TRUE)
  return(this_maf)
}


#' @title Intersect MAF.
#'
#' @description Perform set operations on two MAFs.
#'
#' @details Perform set operations on two MAFs.
#'
#' @param maf1 First list of MAFs.
#' @param maf2 Second list of MAFs.
#' @param set_returned List of MAFs that doesn't share the same start positions as the other list of MAFs. Accepted commands are; "maf1_only" and "maf2_only", default is "maf1_only".
#'
#' @return Set of MAFs with start positions that don't match the start positions in the other supplied MAF file.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' intersected_mafs_l1 = intersect_maf(maf_list1, maf_list2, "maf1_only")
#' intersected_mafs_l2 = intersect_maf(maf_list1, maf_list2, "maf2_only")
#' }
#' 
intersect_maf = function(maf1,
                         maf2,
                         set_returned = "maf1_only"){

  if(set_returned == "maf1_only"){
    maf_set = dplyr::filter(maf1, !Start_Position %in% maf2$Start_Position)
  }else if(set_returned == "maf2_only"){
    maf_set = dplyr::filter(maf2, !Start_Position %in% maf1$Start_Position)
  }
  return(maf_set)
}


#' @title Get Coding SSM Status.
#'
#' @description Tabulate mutation status (SSM) for a set of genes.
#'
#' @details This function takes a vector of gene symbols and subsets the incoming MAF to specified genes. If no genes are provided, the function will default to all lymphoma genes.
#' The function can accept a wide range of incoming MAFs. For example, the user can call this function with `these_samples_metadata` (preferably a metadata table that has been subset to the sample IDs of interest).
#' If this parameter is not called, the function will default to all samples available with [GAMBLR::get_gambl_metadata]. The user can also provide a path to a MAF, or MAF-like file with `maf_path`,
#' or an already loaded MAF can be used with the `maf_data` parameter. If both `maf_path` and `maf_data` is missing, the function will default to run `get_coding_ssm`.
#' This function also has a lot of filtering and convenience parameters giving the user full control of the return. For more information, refer to the parameter descriptions and examples.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR::get_coding_ssm], [GAMBLR::get_ssm_by_patients], [GAMBLR::get_ssm_by_sample],
#' [GAMBLR::get_ssm_by_samples], [GAMBLR::get_ssm_by_region], [GAMBLR::get_ssm_by_regions]
#'
#' @param gene_symbols A vector of gene symbols for which the mutation status will be tabulated. If not provided, lymphoma genes will be returned by default.
#' @param these_samples_metadata The metadata for samples of interest to be included in the returned matrix. Only the column "sample_id" is required. If not provided, the matrix is tabulated for all available samples as default.
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations. Default is TRUE.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param maf_path If the status of coding SSM should be tabulated from a custom maf file, provide path to the maf in this argument. The default is set to NULL.
#' @param maf_data Either a maf loaded from disk or from the database using a get_ssm function.
#' @param include_hotspots Logical parameter indicating whether hotspots object should also be tabulated. Default is TRUE.
#' @param keep_multihit_hotspot Logical parameter indicating whether to keep the gene annotation as mutated when the gene has both hot spot and non-hotspot mutation. Default is FALSE. If set to TRUE, will report the number of non-hotspot mutations instead of tabulating for just mutation presence.
#' @param recurrence_min Integer value indicating minimal recurrence level.
#' @param seq_type The seq_type you want back, default is genome.
#' @param projection Specify projection (grch37 or hg38) of mutations. Default is grch37.
#' @param review_hotspots Logical parameter indicating whether hotspots object should be reviewed to include functionally relevant mutations or rare lymphoma-related genes. Default is TRUE.
#' @param genes_of_interest A vector of genes for hotspot review. Currently only FOXO1, MYD88, and CREBBP are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF file. The default is hg19 genome build.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is TRUE.
#'
#' @return A data frame with tabulated mutation status.
#'
#' @import dplyr tidyr
#' @export
#'
#' @examples
#' coding_tabulated_df = get_coding_ssm_status(maf_data = grande_maf,
#'                                             gene_symbols = "EGFR")
#' 
#' #all lymphoma genes from bundled NHL gene list
#' coding_tabulated_df = get_coding_ssm_status()
#' 
get_coding_ssm_status = function(gene_symbols,
                                 these_samples_metadata,
                                 from_flatfile = TRUE,
                                 augmented = TRUE,
                                 min_read_support = 3,
                                 maf_path = NULL,
                                 maf_data,
                                 include_hotspots = TRUE,
                                 keep_multihit_hotspot = FALSE,
                                 recurrence_min = 5,
                                 seq_type = "genome",
                                 projection = "grch37",
                                 review_hotspots = TRUE,
                                 genes_of_interest = c("FOXO1", "MYD88", "CREBBP"),
                                 genome_build = "hg19",
                                 include_silent = TRUE){

  if(missing(gene_symbols)){
    message("defaulting to all lymphoma genes")
    gene_symbols = pull(lymphoma_genes, Gene)
  }

  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata()
  }

  # call it once so the object can be reused if user wants to annotate hotspots
  if(!missing(maf_data)){
    coding_ssm = maf_data %>%
      dplyr::filter(Variant_Classification %in% coding_class)

  }else if (!is.null(maf_path)){
    coding_ssm = fread_maf(maf_path)
    coding_ssm = coding_ssm %>%
      dplyr::filter(Variant_Classification %in% coding_class)
  }

  if(missing(maf_data) & is.null(maf_path)){
    coding_ssm = get_coding_ssm(projection = projection, seq_type = seq_type, from_flatfile = from_flatfile, augmented = augmented, min_read_support = 3, basic_columns = FALSE, include_silent = include_silent)
  }

  coding_ssm = coding_ssm %>%
    dplyr::filter(Variant_Classification %in% coding_class)

  coding = coding_ssm %>%
    dplyr::filter(Hugo_Symbol %in% gene_symbols & Variant_Classification != "Synonymous") %>%
    dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
    dplyr::rename("sample_id" = "Tumor_Sample_Barcode", "gene" = "Hugo_Symbol") %>%
    unique() %>%
    mutate(mutated = 1)

  samples_table = dplyr::select(these_samples_metadata, sample_id)
  wide_coding = pivot_wider(coding,names_from = "gene", values_from = "mutated",values_fill = 0)
  all_tabulated = left_join(samples_table, wide_coding)
  all_tabulated = all_tabulated %>%
    replace(is.na(.), 0)

  # include hotspots if user chooses to do so
  if(include_hotspots){
    # first annotate
    annotated = annotate_hotspots(coding_ssm, recurrence_min = recurrence_min)
    # review for the supported genes
    if(review_hotspots){
      annotated = review_hotspots(annotated, genes_of_interest = genes_of_interest, genome_build = genome_build)
    }
    message("annotating hotspots")
    hotspots = annotated %>%
      dplyr::filter(Hugo_Symbol %in% genes_of_interest) %>%
      dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, hot_spot) %>%
      dplyr::rename("sample_id" = "Tumor_Sample_Barcode", "gene" = "Hugo_Symbol") %>%
      dplyr::mutate(gene = paste0(gene, "HOTSPOT")) %>%
      unique() %>%
      dplyr::mutate(mutated = ifelse(hot_spot == "TRUE", 1, 0)) %>%
      dplyr::filter(mutated == 1) %>%
      dplyr::select(-hot_spot)

    # long to wide hotspots, samples are tabulated with 0 if no hotspot is detected
    wide_hotspots = pivot_wider(hotspots, names_from = "gene", values_from = "mutated", values_fill = 0)
    # join with the ssm object
    all_tabulated = left_join(all_tabulated, wide_hotspots)
    all_tabulated = all_tabulated %>%
      replace(is.na(.), 0)

    all_tabulated = all_tabulated %>%
      dplyr::select(where(~ any(. != 0)))

    all_tabulated = as.data.frame(all_tabulated)
    # make SSM and hotspots non-redundant by giving priority to hotspot feature and setting SSM to 0
    for (hotspot_site in colnames(wide_hotspots)[grepl("HOTSPOT", colnames(wide_hotspots))]){
      message(hotspot_site)
      this_gene = gsub("HOTSPOT", "", hotspot_site)
      redundant_features = all_tabulated %>%
        dplyr::select(starts_with(this_gene))

      # if not both the gene and the hotspot are present, go to the next iteration
      if(ncol(redundant_features)!= 2) next
      message("OK")
      # if both gene and it's hotspot are in the matrix, give priority to hotspot feature
      all_tabulated[(all_tabulated[, this_gene] >0 & all_tabulated[, paste0(this_gene, "HOTSPOT")] == 1),][,c(this_gene, paste0(this_gene, "HOTSPOT"))][, this_gene] = 0

      # in case gene has both hotspot and another mutation in the same gene,
      # keep both tabulated as multihits
      if(keep_multihit_hotspot){
        # determine which samples have hot spot and another mutation in same gene
        multihits <- annotated %>%
            dplyr::filter(Hugo_Symbol == this_gene) %>%
            group_by(Tumor_Sample_Barcode) %>%
            dplyr::mutate(n_mut = n()) %>%
            dplyr::filter(
                n_mut > 1
            ) %>%
            dplyr::distinct(Tumor_Sample_Barcode, n_mut, hot_spot) %>%
            # account for cases with both hotspot and not hotspot to avoid
            # double-counting the number of mutations
            mutate_at(vars(hot_spot), ~replace_na(., "FALSE")) %>%
            dplyr::mutate(
                n_mut = ifelse(
                    hot_spot == "TRUE",
                    n_mut - 1,
                    n_mut
                )
            ) %>%
            group_by(Tumor_Sample_Barcode) %>%
            dplyr::arrange(n_mut) %>%
            slice_head() %>%
            ungroup %>%
            select(-hot_spot)

        # Return the annotation of this gene to mutated in these samples
        all_tabulated <- all_tabulated %>%
            left_join(
              .,
              multihits,
              by = c("sample_id" = "Tumor_Sample_Barcode")
            ) %>%
            dplyr::mutate(
                {{this_gene}} := ifelse(
                        !is.na(n_mut),
                        n_mut,
                        !!!syms(this_gene)
                    )
            ) %>%
            select(- n_mut)
      }

    }

  }
  return(all_tabulated)
}


#' @title Trim Scale Expressions.
#'
#' @description INTERNAL FUNCTION called by prettyOncoplot, not meant for out-of-package usage.
#'
#' @details INTERNAL FUNCTION called by prettyOncoplot, not meant for out-of-package usage.
#'
#' @param x Numeric value (of expression) to be trimmed.
#'
#' @return Numeric value.
#' 
#' @noRd
#'
#' @examples
#' trimmed = trim_scale_expression(2)
#' 
trim_scale_expression = function(x){
  quants = unname(quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))
  x = ifelse(x < quants[1], quants[1], x)
  x = ifelse(x > quants[2], quants[2], x)
  x = (x - quants[1]) / (quants[2] - quants[1])
  return(x)
}



#' @title Calculate Mutation Frequency By Sliding Window.
#'
#' @description Count the number of mutations in a sliding window across a region for all samples.
#'
#' @details This function is called to return the mutation frequency for a given region, for all GAMBL samples. Regions are specified with the `this_region`parameter.
#' Alternatively, the region of interest can also be specified by calling the function with `chromosome`, `start_pos`, and `end_pos` parameters.
#' It is also possible to return a plot of the created bins. This is done with setting `plot_type = TRUE`.
#' There are a collection of parameters available for further customizing the return, for more information, refer to the parameter descriptions and examples.
#' This function is unlikely to be used directly in most cases. See [GAMBLR::get_mutation_frequency_bin_matrix] instead.
#'
#' @param this_region Genomic region in bed format.
#' @param chromosome Chromosome name in region.
#' @param start_pos Start coordinate of region.
#' @param end_pos End coordinate of region.
#' @param metadata Data frame containing sample ids and column with annotated data for the 2 groups of interest. All other columns are ignored. Currently, function exits if asked to compare more than 2 groups.
#' @param seq_type The seq_type you want back, default is genome.
#' @param slide_by Slide size for sliding window, default is 100.
#' @param window_size Size of sliding window, default is 1000.
#' @param plot_type Set to TRUE for a plot of your bins. By default no plots are made.
#' @param sortByColumns Which of the metadata to sort on for the heatmap
#' @param return_format Return format of mutations. Accepted inputs are "long" and "long-simple". Default is "long-simple".
#' @param min_count_per_bin Minimum counts per bin, default is 3.
#' @param return_count Boolean statement to return count. Default is FALSE.
#' @param drop_unmutated This may not currently work properly. Default is FALSE.
#' @param classification_column Only used for plotting, default is "lymphgen"
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flat-files (only works for streamlined data, not full MAF details). Default is FALSE.
#' @param mode Only works with indexed flat-files. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#' @return Count matrix.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr tidyr cowplot ggplot2
#' @export
#'
#' @examples
#' chr11_mut_freq = calc_mutation_frequency_sliding_windows(this_region = "chr11:69455000-69459900",
#'                                                          slide_by = 10,
#'                                                          window_size = 10000)
#'
calc_mutation_frequency_sliding_windows = function(this_region,
                                                   chromosome,
                                                   start_pos,
                                                   end_pos,
                                                   metadata,
                                                   seq_type,
                                                   slide_by = 100,
                                                   window_size = 1000,
                                                   plot_type = "none",
                                                   sortByColumns = "pathology",
                                                   return_format = "long-simple",
                                                   min_count_per_bin = 3,
                                                   return_count = FALSE,
                                                   drop_unmutated = FALSE,
                                                   classification_column = "lymphgen",
                                                   from_indexed_flatfile = FALSE,
                                                   mode = "slms-3"){

  max_region = 1000000
  if(missing(metadata)){
    metadata = get_gambl_metadata()
  }
  if(missing(this_region)){
    this_region = paste0(chromosome, ":", start_pos, "-", end_pos)
  }else{
    chunks = region_to_chunks(this_region)
    chromosome = chunks$chromosome
    start_pos = as.numeric(chunks$start)
    end_pos = as.numeric(chunks$end)
  }
  region_size = end_pos - start_pos
  if(region_size < max_region){
    message(paste("processing bins of size", window_size, "across", region_size, "bp region"))
  }else{
    message(paste("CAUTION!\n", region_size, "exceeds maximum size recommended by this function."))
  }
  windows = data.frame(start = seq(start_pos, end_pos, by = slide_by)) %>%
    mutate(end = start + window_size - 1)
  #use foverlaps to assign mutations to bins
  windows.dt = as.data.table(windows)
  region_ssm = GAMBLR::get_ssm_by_region(region = this_region, streamlined = FALSE, seq_type=seq_type, from_indexed_flatfile = from_indexed_flatfile, mode = mode) %>%
    dplyr::rename(c("start" = "Start_Position", "sample_id" = "Tumor_Sample_Barcode")) %>%
    mutate(mutated = 1)

  region.dt = region_ssm %>%
    dplyr::mutate(start = as.numeric(as.character(start)), end = start + 1, end = as.numeric(as.character(end))) %>%
    dplyr::relocate(start, .before=end) %>%
    as.data.table()

  setkey(windows.dt, start, end)
  setkey(region.dt, start, end)
  windows_overlap = foverlaps(windows.dt, region.dt) %>%
    dplyr::filter(!is.na(start)) %>%
    dplyr::rename(c("window_start" = "i.start", "mutation_position" = "start")) %>%
    dplyr::select(-i.end, -end, -mutation_position) %>%
    as.data.frame()

  windows_tallied_full = windows_overlap %>%
    group_by(sample_id, window_start) %>%
    tally() %>%
    dplyr::filter(n >= min_count_per_bin) %>%
    arrange(sample_id) %>%
    as.data.frame()
  windows_tallied = windows_tallied_full

  all_samples = pull(metadata, sample_id) %>%
    unique()

  num_samples = length(all_samples)
  lg_cols = get_gambl_colours("lymphgen")
  path_cols = get_gambl_colours("pathology")
  annos = data.frame(window_start = rep(start_pos,num_samples), sample_id = factor(all_samples))
  annos = left_join(annos, metadata, by = "sample_id")
  windows_tallied = left_join(metadata, windows_tallied, by = "sample_id")
  windows_tallied = arrange(windows_tallied, lymphgen)
  windows_tallied$classification = factor(windows_tallied[,classification_column], levels = unique(windows_tallied[,classification_column]))
  if(drop_unmutated){
    windows_tallied = windows_tallied %>%
      dplyr::filter(!is.na(n))
  }
  if(classification_column == "lymphgen"){
    windows_tallied = arrange(windows_tallied, pathology, lymphgen)
    annos = arrange(annos, pathology, lymphgen)
  }else{
    windows_tallied = arrange(windows_tallied, classification)
    annos = arrange(annos, classification)
  }
  annos$sample_id = factor(annos$sample_id, levels = unique(annos$sample_id))
  windows_tallied$sample_id = factor(windows_tallied$sample_id, levels = unique(windows_tallied$sample_id))
  if(plot_type %in% c("points", "point")){
    #add a bin at position 1 for pathology
    windows_tallied = dplyr::filter(windows_tallied, !is.na(window_start))
    p = ggplot2::ggplot() +
                 geom_point(data = annos, aes(x = window_start, y = sample_id, colour = pathology)) +
                 geom_point(data = windows_tallied, aes(x = window_start, y = sample_id, colour = classification)) +
                 theme(axis.text = element_text(size = 4)) +
                 theme(axis.text.y = element_blank()) +
                 scale_colour_manual(values = c(lg_cols, path_cols))
  }else if(plot_type %in% c("tile", "tiled")){
    print(annos)
    p = ggplot() +
        geom_point(data = annos, aes(x = window_start, y = sample_id, colour = lymphgen)) +
        geom_tile(data = windows_tallied, aes(x = window_start, y = sample_id, fill = n)) +
        scale_fill_gradient(low = "orange", high = "red", na.value = NA) +
        theme_cowplot() +
        scale_colour_manual(values = c(lg_cols, path_cols)) +
        theme(axis.text.y = element_blank())
  }
  if(plot_type != "none"){
    print(p)
  }
  if(return_count){
    windows_tallied = mutate(windows_tallied, bin = paste0(window_start, "_", chromosome)) %>%
      mutate(mutated = n)
  }else{
    #binary mutated or not
    windows_tallied = mutate(windows_tallied, bin = paste0(window_start, "_", chromosome)) %>%
    mutate(mutated = 1)
  }
  a = dplyr::select(windows_tallied, sample_id, bin, mutated)
  completed = complete(a, sample_id, bin, fill = list(mutated = 0))
  widened = pivot_wider(completed, names_from = sample_id, values_from = mutated)
  if(return_format == "long"){
    return(windows_tallied)
  }else if(return_format == "long-simple"){
    #just return the columns needed for completing and making a wide matrix for many regions
    windows_simple = dplyr::select(windows_tallied, sample_id, bin, mutated)
    return(windows_simple)
  }else{
    return(widened)
  }
}


#' @title SV To BED File.
#'
#' @description Write bedpe format data frame to a file that will work with IGV and UCSC genome browser.
#'
#' @details This function takes four parameters; a data frame with SVs, formatted as bedpe with `sv_df`.
#' The `path` parameter lets the user control the output folder. The default is "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/misc/".
#' `file_name` for specifying the output file name. Results will be written to `results/icgc_dart/misc/`.
#' Lastly, the `add_chr_prefix` lets the user control if chromosomes should be prefixed with "chr" or not.
#' The default is TRUE.
#'
#' @param sv_df data frame of bedpe formatted SV data.
#' @param path The path to the output folder. Default is "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/misc/".
#' @param filename File name (will be written to results/icgc_dart/misc/FILENAME).
#' @param add_chr_prefix Whether to force chr to be added to chromosome names. Default is TRUE.
#'
#' @return bedpe data frame that is compatible with IGN browser.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' SVs_bedpe = sv_to_bedpe_file(sv_df = sv_dataframe, 
#'                              filename = "SVs.bedpe",
#'                              add_chr_prefix = TRUE)
#' }
#' 
sv_to_bedpe_file = function(sv_df,
                            path = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/misc/",
                            filename = "my_svs.bedpe",
                            add_chr_prefix = TRUE){

  #add chr prefix if missing
  if(add_chr_prefix){
    if(!any(grepl("chr", region_sv$CHROM_A[1]))){
      sv_df = mutate(sv_df, CHROM_A = paste0("chr", CHROM_A)) %>%
        mutate(CHROM_B = paste0("chr", CHROM_B))
    }
  }
  bed_file = paste0(path, filename)
  write.table(sv_df, file = bed_file, sep = "\t", quote = F, row.names = F, col.names = F)
}


#' @title Region To Chunks.
#'
#' @description Parse a region string into; chromosome, start and end.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::calc_mutation_frequency_sliding_windows], not meant for out-of-package usage.
#'
#' @param region A region string e.g. "chrX:12345-678910".
#'
#' @return A named list.
#' 
#' @noRd
#'
#' @examples
#' chr_start_end = region_to_chunks("chr1:1111-2222")
#' 
region_to_chunks = function(region){

  region = unname(region)
  region = gsub(",", "", region)
  #format is chr6:37060224-37151701
  split_chunks = unlist(strsplit(region, ":"))
  chromosome = split_chunks[1]
  startend = unlist(strsplit(split_chunks[2], "-"))
  qstart = startend[1]
  qend = startend[2]
  return(list(chromosome = chromosome, start = qstart, end = qend))
}


#' @title Convert mutation data to a shareable format.
#'
#' @description `sanitize_maf_data` returns an oncomatrix of patient/gene data indicating only data needed to produce an oncoplot.
#'
#' @details Write an oncomatrix from a MAF File for further plotting. This is meant to be run by individuals who have access to data sets to
#' "sanitize" a subset of data for subsequent use by them or others who don't have permission to access the raw data.
#' Example: User J has full permissions for ICGC data and has read permissions on a MAF file. User B needs to make some oncoplots
#' and/or perform some statistical analysis on the frequency and assortment of mutations in that data set but don't need all the details.
#' User J can run this function on a maf file and provide the path of the output to user B.
#'
#' @param mutation_maf_path Provide either the full path to a MAF file.
#' @param mutation_maf_data Otherwise provide a data frame of the MAF data.
#' @param output_oncomatrix Optionally provide the path for your sanitized output file (otherwise it writes to the working directory).
#' @param genes_keep Specify which genes you want to remain in the output. Make sure there are no duplicated elements in the vector.
#' @param genes_drop Optionally specify which genes to drop (this doesn't mean all other genes will remain. Maftools decides that part).
#'
#' @return The full path to the oncomatrix file (a matrix with Variant_Classification or Multi_Hit indicating coding mutation status per patient).
#'
#' @import maftools dplyr
#' @export
#'
#' @examples
#' 
#' safe_oncomatrix_path = sanitize_maf_data(mutation_maf_data = grande_maf, 
#'                                          genes_keep = c("MYC", "ID3", "ARID1A",
#'                                                         "FOXO1", "TP53", "FAT4",
#'                                                         "IGLL5"))
#'
sanitize_maf_data = function(mutation_maf_path,
                             mutation_maf_data,
                             output_oncomatrix,
                             genes_keep,
                             genes_drop = c()){

  if(missing(mutation_maf_path) & missing(mutation_maf_data)){

    warning("Provide either a path to a MAF file or a data frame of mutations")
    return()
  }
  if(!missing(mutation_maf_path)){
    mutation_maf_data = fread_maf(mutation_maf_path)
  }
  #because we use oncoplot we need to ensure the user gets an oncomatrix for all the genes they care about
  #optionally we also exclude genes
  if(missing(genes_keep)){
    warning("you should provide a list of genes to retain in the output to ensure your genes of interest are included")
    genes_keep = rev(unique(lymphoma_genes$Gene))
  }
  maf_o = maftools::read.maf(mutation_maf_data)

  maftools::oncoplot(maf_o, genes = genes_keep, writeMatrix = T, removeNonMutated = F)  #writes to working directory
  if(!missing(output_oncomatrix)){
    #rename it
    file.rename("onco_matrix.txt", output_oncomatrix)
  }else{
    output_oncomatrix = paste0(getwd(), "/onco_matrix.tsv")
  }
  message(paste("your data is in:", output_oncomatrix))
  return(output_oncomatrix)
}


#' @title Annotate Hotspots.
#'
#' @description Annotate MAF-like data frome with a hot_spot column indicating recurrent mutations.
#'
#' @details This function takes an already loaded MAF data frame with the `mutation_maf` parameter.
#' The user can then control the minimum number of recurrences for mutations to be included with `recurrance_min`,
#' The default is 5. `analysis_base` controls the base name go hotspot output directory.
#' Lastly, `p_thresh` sets the p value threshold, default is 0.05.
#'
#' @param mutation_maf A data frame in MAF format.
#' @param recurrence_min minimum number of recurrences for mutation to be included, default is 5.
#' @param analysis_base Base name for hot spot output directory.
#' @param p_thresh P value threshold, default is 0.05.
#'
#' @return The same data frame with one additional column "hot_spot".
#'
#' @import dplyr tidyr readr
#' @export
#'
#' @examples
#' my_metadata = get_gambl_metadata()
#' all_coding_ssm = get_coding_ssm(these_samples_metadata = my_metadata,
#'                                 projection = "grch37", 
#'                                 seq_type = "genome")
#'
#' hot_ssms = annotate_hotspots(all_coding_ssm)
#'
annotate_hotspots = function(mutation_maf,
                             recurrence_min = 5,
                             analysis_base = c("FL--DLBCL", "BL--DLBCL"),
                             p_thresh = 0.05){

  hotspot_info = list()
  for(abase in analysis_base){
    base_path = check_config_value(config::get("repo_base"))

    clust_full_path = paste0(base_path, check_config_value(config::get("results_versioned")$oncodriveclustl$clusters))
    clust_full_path = glue::glue(clust_full_path)
    all_full_path = paste0(base_path, check_config_value(config::get("results_versioned")$oncodriveclustl$elements))
    all_full_path = glue::glue(all_full_path)
    clust_hotspot = suppressMessages(readr::read_tsv(clust_full_path))
    all_hotspot = suppressMessages(readr::read_tsv(all_full_path))

  clustered_hotspots = clust_hotspot %>%
    dplyr::select(-RANK) %>%
    dplyr::filter(N_SAMPLES > recurrence_min & P < p_thresh)

    arranged = clustered_hotspots %>%
      separate_rows(COORDINATES, convert = TRUE) %>%
      group_by(SYMBOL, MAX_COORD) %>%
      arrange(COORDINATES)

    mins = arranged %>%
      slice_head() %>%
      dplyr::rename("START" = "COORDINATES")

    maxs = arranged %>%
      slice_tail() %>%
      dplyr::rename("END" = "COORDINATES")

    hotspot_ranges = left_join(mins, dplyr::select(maxs, c(MAX_COORD, END)), by = c("SYMBOL", "MAX_COORD"))
    hotspot_info[[abase]] = hotspot_ranges
  }
  merged_hotspot = do.call("rbind", hotspot_info) %>%
    ungroup()

  long_hotspot = merged_hotspot %>%
    dplyr::select(MAX_COORD, CHROMOSOME, START, END) %>%
    pivot_longer(c(START, END), names_to = "which", values_to = "COORDINATE") %>%
      dplyr::select(-which)

  #again take highest and lowest value for each MAX_COORD
  starts = long_hotspot %>%
    group_by(MAX_COORD) %>%
    arrange(COORDINATE) %>%
    slice_head()

  ends = long_hotspot %>%
    group_by(MAX_COORD) %>%
    arrange(COORDINATE) %>%
    slice_tail()

  long_hotspot = bind_rows(starts, ends)
  filled_coords = long_hotspot %>%
    group_by(MAX_COORD) %>%
    arrange(MAX_COORD, COORDINATE) %>%
    complete(COORDINATE = seq(COORDINATE[1], COORDINATE[2])) %>%
    tidyr::fill(CHROMOSOME, .direction = "up") %>%
    dplyr::rename("Start_Position" = "COORDINATE") %>%
    dplyr::rename("Chromosome" = "CHROMOSOME") %>%
    ungroup()

  filled_coords = mutate(filled_coords, hot_spot = TRUE)
  #just the ssms that match these coordinates!
  hot_ssms = left_join(mutation_maf, filled_coords, by = c("Chromosome", "Start_Position"))
  return(hot_ssms)
}


#' @title Review Hotspots.
#'
#' @description Annotate MAF-like data frome with a hot_spot column indicating recurrent mutations.
#'
#' @details This function takes an annotated MAF (with [GAMBLR::annotate_hotspots]) and adds a new column, "hot_spot", to the same data frame.
#' Genes for hotspot review are supplied with the `genes_of_interest` parameter.
#' Currently only a few sets of genes are supported, see parameter description for more information and limitations.
#' The desired genome build can be specified with `genome_build` parameter. Should be the same as the incoming MAF.
#'
#' @param annotated_maf A data frame in MAF format that has hotspots annotated using the function annotate_hotspots().
#' @param genes_of_interest A vector of genes for hotspot review. Currently only FOXO1, MYD88, CREBBP, NOTCH1, NOTCH2, CD79B and EZH2 are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF file. The default is grch37 genome build.
#'
#' @return The same data frame (as given to the `annotated_maf` parameter) with the reviewed column "hot_spot".
#'
#' @import dplyr
#' @export
#'
#' @examples
#' hot_ssms = review_hotspots(annotate_hotspots(get_coding_ssm(seq_type = "genome")),
#'                            genes_of_interest = c("CREBBP"))
#'
review_hotspots = function(annotated_maf,
                           genes_of_interest = c("FOXO1", "MYD88", "CREBBP", "NOTCH1", "NOTCH2", "CD79B", "EZH2"),
                           genome_build = "grch37"){

  # define the list of genes currently supported for review
  supported_genes = c("FOXO1", "MYD88", "CREBBP", "NOTCH1", "NOTCH2", "CD79B", "EZH2")

  # check genome build because CREBBP coordinates are hg19-based or hg38-based

  if (genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
    coordinates = GAMBLR::hotspot_regions_grch37
  }else if(genome_build %in% c("hg38", "grch38", "GRCh38")){
    coordinates = GAMBLR::hotspot_regions_hg38
  }else{
    stop("The genome build specified is not currently supported. Please provide MAF file in one of the following cordinates: hg19, grch37, hs37d5, GRCh37, hg38, grch38, or GRCh38")
  }
  # check that at least one of the currently supported genes are present
  if (length(intersect(supported_genes, genes_of_interest))==0){
      stop(paste0("Currently only ",  paste(supported_genes, collapse=", "), " are supported. Please specify one of these genes."))
  }
  # notify user that there is limited number of genes currently supported
  if (length(setdiff(genes_of_interest, supported_genes))>0){
      message(strwrap(paste0("Currently only ", paste(supported_genes, collapse=", "),
                             " are supported. By default only these genes from the supplied list will be reviewed. Reviewing hotspots for genes ",
                             paste(intersect(supported_genes, genes_of_interest), collapse = ", "), ", it will take a second ...")))
  }
  if("FOXO1" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "FOXO1" &
                                        HGVSp_Short == "p.M1?",
                                        "TRUE", hot_spot))
  }

  if("CREBBP" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CREBBP" &
                                        Start_Position > coordinates["CREBBP", "start"] &
                                        End_Position < coordinates["CREBBP", "end"] &
                                        Variant_Classification == "Missense_Mutation",
                                        "TRUE", hot_spot))
  }
  if("EZH2" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "EZH2" &
                                        Start_Position > coordinates["EZH2", "start"] &
                                        End_Position < coordinates["EZH2", "end"],
                                        "TRUE", hot_spot))
  }
  if("MYD88" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "MYD88" &
                                        HGVSp_Short %in% c("p.L273P", "p.L265P"),
                                        "TRUE", hot_spot))
  }
  if("NOTCH1" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "NOTCH1" &
                                        Start_Position < coordinates["NOTCH1", "start"],
                                        "TRUE", hot_spot))
  }
  if("NOTCH2" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "NOTCH2" &
                                        Start_Position < coordinates["NOTCH2", "start"],
                                        "TRUE", hot_spot))
  }

  if("CD79B" %in% genes_of_interest){
      truncating_variants = c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Region", "Splice_Site")
      annotated_maf = annotated_maf %>%
         dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CD79B" &
                                         Start_Position < coordinates["CD79B_trunc", "start"] &
                                         Variant_Classification %in% truncating_variants,
                                         "TRUE", hot_spot)) %>%
          dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CD79B" &
                                          Start_Position < coordinates["CD79B_NONtrunc", "start"] &
                                          ! Variant_Classification %in% truncating_variants,
                                          "TRUE", hot_spot))
  }
  return(annotated_maf)
}


#' @title SV To Custom Track.
#'
#' @description Make a UCSC-ready custom track file from SV data.
#'
#' @details This function takes an incoming SV data frame and outputs a bed file, ready for visualization on the UCSC Genome Browser.
#' Specify the output file with `output_file`, indicate if the incoming SVs are annotated with `is_annotated` (default is TRUE).
#' Lastly, the user can also specify if the incoming SV data frame should be subset to specific mutation types (e.g deletions, duplications, insertions, etc.).
#' This is specified with the `sv_name` parameter. Default is to include all SV subtypes.
#'
#' @param sv_bedpe A bedpe formatted data frame of SVs.
#' @param output_file A bed file with UCSC custom header.
#' @param is_annotated Set to TRUE if input SV bedpe is annotated, default is TRUE.
#' @param sv_name SV name. Default is set to "all" = include all subtypes of SVs.
#'
#' @return Nothing.
#'
#' @import dplyr tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' #custom track with annotations
#' all_sv = get_manta_sv(verbose = FALSE)
#' annotated_sv = annotate_sv(sv_data = all_sv)
#' sv_to_custom_track(annotated_sv,
#'                    output_file = "GAMBL_sv_custom_track_annotated.bed",
#'                    is_annotated = TRUE)
#'
#' #custom track (no anotatated SVs)
#' sv_to_custom_track(all_sv,
#'                    output_file = "GAMBL_sv_custom_track_annotated.bed",
#'                    is_annotated = FALSE)
#' }
#' 
sv_to_custom_track = function(sv_bedpe,
                              output_file,
                              is_annotated = TRUE,
                              sv_name = "all"){

  if(is_annotated){
    #reduce to a bed-like format
    sv_data1 = mutate(annotated_sv, annotation = paste0(chrom1, ":", start1, "_", fusion)) %>%
      dplyr::select(chrom2, start2, end2, tumour_sample_id, annotation, fusion)

    sv_data2 = mutate(annotated_sv, annotation = paste0(chrom2, ":", start2, "_", fusion)) %>%
      dplyr::select(chrom1, start1, end1, tumour_sample_id, annotation, fusion)

    print(head(sv_data1))
    print(head(sv_data2))
    colnames(sv_data1) = c("chrom", "start", "end", "sample_id", "annotation", "fusion")
    colnames(sv_data2) = c("chrom", "start", "end", "sample_id", "annotation", "fusion")
    sv_data = bind_rows(sv_data1, sv_data2)
    sv_data = mutate(sv_data, end = end + 10)
  }else{
    sv_data_1 = mutate(sv_bedpe, annotation = paste0( CHROM_B, ":", START_B)) %>%
      dplyr::select(CHROM_A, START_A, END_A, tumour_sample_id, annotation)

    sv_data_2 = mutate(sv_bedpe, annotation = paste0( CHROM_A, ":", START_A)) %>%
      dplyr::select(CHROM_B, START_B, END_B, tumour_sample_id, annotation)

    colnames(sv_data_1) = c("chrom", "start", "end", "sample_id", "annotation")
    colnames(sv_data_2) = c("chrom", "start", "end", "sample_id", "annotation")
    sv_data = bind_rows(sv_data_1, sv_data_2)

  }
  if(!any(grepl("chr", sv_data[,1]))){
    #add chr
    sv_data[,1] = unlist(lapply(sv_data[,1], function(x){paste0("chr", x)}))
  }

  coo_cols = get_gambl_colours("COO")
  path_cols = get_gambl_colours("pathology")
  all_cols = c(coo_cols, path_cols)
  colour_df = data.frame(coo = names(all_cols), colour = all_cols)

  rgb_df = data.frame(t(col2rgb(all_cols))) %>%
    mutate(consensus_coo_dhitsig = names(all_cols)) %>%
    unite(col = "rgb", red, green, blue, sep = ",")

  meta = get_gambl_metadata() %>%
    dplyr::select(sample_id, "consensus_coo_dhitsig", pathology) %>%
    mutate(consensus_coo_dhitsig = if_else(consensus_coo_dhitsig == "NA", pathology, consensus_coo_dhitsig))

  samples_coloured = left_join(meta, rgb_df)
  sv_bed_coloured = left_join(sv_data, samples_coloured) %>%
    arrange(pathology)

  write_bed = function(coloured_svs, sv_name, output_file_base){
    data_bed = coloured_svs %>%
      mutate(details = paste0(sample_id, "_", annotation)) %>%
      mutate(score = 0, strand = "+", end = end + 1, start1 = start, end1 = end) %>%
      dplyr::select(chrom, start, end, details, score, strand, start1, end1, rgb) %>%
      dplyr::filter(!is.na(rgb)) %>%
      unique()

    header_content = paste0('track name="GAMBL SVs ', sv_name, '" description="SV breakpoints ', sv_name, '" visibility=2 itemRgb="On"\n')
    cat(header_content, file = output_file)
    message(paste("writing to", output_file))
    tabular = write.table(data_bed, file = output_file, quote = F, sep = "\t", row.names = F, col.names = F, append = TRUE)
  }
  write_bed(sv_bed_coloured, sv_name = sv_name)
}


#' @title Maf To Custom Track.
#'
#' @description Convert a maf-formatted data frame into a bed custom track file for UCSC.
#'
#' @details This function takes an incoming MAF and converts it to a UCSC Genome Browser ready BED (or bigbed/biglolly) file.
#' Optional parameters available for further customization of the returned file. For more information, refer to the parameter descriptions and function examples.
#'
#' @param maf_data Either a maf loaded from disk or from the database using a get_ssm function.
#' @param these_samples_metadata Optional argument, a metadata table subset to the samples of interest. If not provided, the function will run [GAMBLR::get_gambl_metadata] for all available samples.
#' @param seq_type The seq type you want back, default is "genome".
#' @param output_file Name for your new bed file that can be uploaded as a custom track to UCSC.
#' @param as_bigbed Boolean parameter controlling the format of the returned file. Default is FALSE.
#' @param colour_column Set the colouring properties of the returned bed file. Per default, this function will assign colour based on "lymphgen".
#' @param as_biglolly Boolean parameter controlling the format of the returned file. Default is FALSE (i.e a BED file will be returned).
#' @param track_name Track name. Default is "GAMBL mutations"
#' @param track_description Track description. Default is "mutations from GAMBL"
#' @param verbose Default is FALSE.
#' @param padding_size Optional parameter specifying the padding size in the returned file, default is 0.
#'
#' @return Nothing.
#'
#' @import tidyr dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' maf_to_custom_track(my_maf_data, "/home/rmorin/private/some_mutations.bed")
#' }
#' 
maf_to_custom_track = function(maf_data,
                               these_samples_metadata,
                               seq_type = "genome",
                               output_file,
                               as_bigbed = FALSE,
                               colour_column = "lymphgen",
                               as_biglolly = FALSE,
                               track_name = "GAMBL mutations",
                               track_description = "mutations from GAMBL",
                               verbose = FALSE,
                               padding_size = 0){

  #reduce to a bed-like format
  maf_data = dplyr::select(maf_data, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)
  colnames(maf_data) = c("chrom", "start", "end", "sample_id")
  maf_data = mutate(maf_data,end = end + padding_size)
  if(!any(grepl("chr", maf_data[,1]))){
    #add chr
    maf_data[,1] = unlist(lapply(maf_data[,1], function(x){paste0("chr", x)}))
  }
  lymphgen_cols = get_gambl_colours(colour_column,verbose=verbose)

  colour_df = data.frame(group = names(lymphgen_cols), colour = lymphgen_cols)

  rgb_df = data.frame(t(col2rgb(lymphgen_cols))) %>%
    mutate(group = names(lymphgen_cols),hex=unname(lymphgen_cols)) %>%
    unite(col = "rgb", red, green, blue, sep = ",")
  if(verbose){
    print(rgb_df)
  }
  if(missing(these_samples_metadata)){
    meta = get_gambl_metadata(seq_type_filter = seq_type) %>% dplyr::select(sample_id,all_of(colour_column))
  }else{
    meta = these_samples_metadata %>% dplyr::select(sample_id,all_of(colour_column))
  }
  colnames(meta)[2]="group"


  samples_coloured = left_join(meta, rgb_df)
  if(verbose){
    print(samples_coloured)
  }

  maf_bed = maf_data %>%
    mutate(score = 0, strand = "+", start1 = start-1,start=start1, end1 = end)
  if(verbose){
    print(head(maf_bed))
  }
  maf_coloured = left_join(maf_bed, samples_coloured, by = "sample_id") %>%
    dplyr::select(-group) %>%
    mutate(rgb=ifelse(is.na(rgb),"0,0,0",rgb))
  maf_summary = group_by(maf_coloured,hex) %>% tally()
  if(verbose){
    print(maf_summary)
    print(head(maf_coloured))
  }
  maf_coloured = dplyr::select(maf_coloured,-hex)
  if(as_bigbed | as_biglolly){

    if(grepl(pattern = ".bb$",x = output_file)){
      #temp file will be .bed
      temp_bed = gsub(".bb$",".bed",output_file)

    }else{
      stop("please provide an output file name ending in .bb to create a bigBed file")
    }

    maf_coloured = mutate(maf_coloured,sample_id="redacted") %>%
      arrange(chrom,start)
    if(as_biglolly){
      #currently the same code is run either way but this may change so I've separated this until we settle on format
      #TO DO: collapse based on hot spot definition and update column 4 (score) based on recurrence
      #needs to have size column
      maf_score_options = factor(maf_coloured$rgb)
      maf_coloured$score = as.numeric(maf_score_options)

      #determine frequency of each event per group to assign the size
      maf_coloured = group_by(maf_coloured,start,rgb) %>% mutate(size=n())

      #maf_coloured = mutate(maf_coloured,size=10)

      write.table(maf_coloured, file = temp_bed, quote = F, sep = "\t", row.names = F, col.names = F)
      #conversion:
      autosql_file = "/Users/rmorin/git/LLMPP/resources/reference/ucsc/bigLollyExample3.as"

      bigbedtobed = "/Users/rmorin/miniconda3/envs/ucsc/bin/bedToBigBed"
      bigbed_conversion = paste0(bigbedtobed," -as=",autosql_file," -type=bed9+1 ",temp_bed," /Users/rmorin/git/LLMPP/resources/reference/ucsc/hg19.chrom.sizes ",output_file)
      print(bigbed_conversion)
      system(bigbed_conversion)
    }else{
      write.table(maf_coloured, file = temp_bed, quote = F, sep = "\t", row.names = F, col.names = F)
      #conversion:
      bigbedtobed = "/Users/rmorin/miniconda3/envs/ucsc/bin/bedToBigBed"
      bigbed_conversion = paste(bigbedtobed,"-type=bed9",temp_bed,"/Users/rmorin/git/LLMPP/resources/reference/ucsc/hg19.chrom.sizes",output_file)

      system(bigbed_conversion)
    }
  }else{
    header_ucsc = paste0('track name="',track_name,'" description="', track_description, '" visibility=2 itemRgb="On"\n')
    cat(header_ucsc,file = output_file)
    write.table(maf_coloured, file = output_file, quote = F, sep = "\t", row.names = F, col.names = F, append = TRUE)
  }
}

test_glue = function(placeholder="INSERTED"){
  some_string = "this text has {placeholder}"
  print(some_string)
  ss=glue::glue(some_string)
  print(ss)
}


#' @title Collate Results
#'
#' @description Bring together all derived sample-level results from many GAMBL pipelines.
#'
#' @details This function takes a data frame with sample IDs (in the first column) with the `sample_table` parameter and adds sample-level results from many of the available GAMBL pipelines.
#' Optional parameters are `these_samples_metadata` and `join_with_full_metadata`. If `join_with_full_metadata` is set to TRUE, the function can either work with an already subset metadata
#' table (`these_sampels_metadata`), or, if not provided, the function will default to all metadata returned with `get_gambl_metadata`, allowing the user to extend the available information in a metadata table.
#' This function has also been designed so that it can get cached results, meaning that not all individual collate helper functions would have to be run to get results back.
#' To do so, run this function with `from_cache = TRUE` (default). In addition, it's also possible to regenerate the cached results, this is done by setting `write_to_file = TRUE`,
#' This operation auto defaults `from_cache = FALSE`. `case_set` is an optional parameter available for subsetting the return to an already defined set of cases.
#' Lastly, `seq_type_filter` lets the user control what seq type results will be returned for. Default is "genome". For more information on how to get the most out of this function,
#' refer to function examples, vignettes and parameter descriptions.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param write_to_file Boolean statement that outputs tsv file (/projects/nhl_meta_analysis_scratch/gambl/results_local/shared/gambl_{seq_type_filter}_results.tsv) if TRUE, default is FALSE.
#' @param join_with_full_metadata Join with all columns of metadata, default is FALSE.
#' @param these_samples_metadata Optional argument to use a user specified metadata df, overwrites get_gambl_metadata in join_with_full_metadata.
#' @param case_set Optional short name for a pre-defined set of cases.
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param from_cache Boolean variable for using cached results (/projects/nhl_meta_analysis_scratch/gambl/results_local/shared/gambl_{seq_type_filter}_results.tsv), default is TRUE. If write_to_file is TRUE, this parameter auto-defaults to FALSE.
#'
#' @return A table keyed on biopsy_id that contains a bunch of per-sample results from GAMBL
#'
#' @import dplyr readr config
#' @export
#'
#' @examples
#' library(dplyr)
#' 
#' #get collated results for all capture samples, using cached results
#' capture_collated_everything = collate_results(seq_type_filter = "capture",
#'                                               from_cache = TRUE,
#'                                               write_to_file = FALSE)
#'
#' #use an already subset metadata table for getting collated results (cached)
#' metadata = get_gambl_metadata(seq_type_filter = "genome") %>%
#'  dplyr::filter(pathology == "FL")
#' 
#' fl_collated = collate_results(seq_type_filter = "genome",
#'                               join_with_full_metadata = TRUE,
#'                               these_samples_metadata = metadata,
#'                               write_to_file = FALSE,
#'                               from_cache = TRUE)
#'
#' #get collated results for all genome samples and join with full metadata
#' everything_collated = collate_results(seq_type_filter = "genome",
#'                                       from_cache = TRUE,
#'                                       join_with_full_metadata = TRUE)
#'
#' #another example demonstrating correct usage of the sample_table parameter.
#' fl_samples = get_gambl_metadata(seq_type_filter = "genome") %>%
#'  dplyr::filter(pathology == "FL") %>%
#'  dplyr::select(sample_id, patient_id, biopsy_id)
#' 
#' fl_collated = collate_results(sample_table = fl_samples,
#'                               seq_type_filter = "genome",
#'                               from_cache = TRUE)
#'
collate_results = function(sample_table,
                           write_to_file = FALSE,
                           join_with_full_metadata = FALSE,
                           these_samples_metadata,
                           case_set,
                           sbs_manipulation = "",
                           seq_type_filter = "genome",
                           from_cache = TRUE){

  # important: if you are collating results from anything but WGS (e.g RNA-seq libraries) be sure to use biopsy ID as the key in your join
  # the sample_id should probably not even be in this file if we want this to be biopsy-centric
  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter = seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  if(write_to_file){
    from_cache = FALSE #override default automatically for nonsense combination of options
  }

  #get paths to cached results, for from_cache = TRUE and for writing new cached results.
  output_file = check_config_value(config::get("results_merged")$collated)
  output_base = check_config_value(config::get("project_base"))
  output_file = paste0(output_base, output_file)
  output_file = glue::glue(output_file)
  print(output_file)
  if(from_cache){
    #check for missingness
    if(!file.exists(output_file)){
      print(paste("missing: ", output_file))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    #read cached results
    sample_table = suppressMessages(read_tsv(output_file) %>% dplyr::filter(sample_id %in% sample_table$sample_id))

  }else{
    message("Slow option: not using cached result. I suggest from_cache = TRUE whenever possible")
    #edit this function and add a new function to load any additional results into the main summary table
    sample_table = collate_ssm_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_sv_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_curated_sv_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_ashm_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_nfkbiz_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_csr_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_ancestry(sample_table = sample_table, seq_type_filter = seq_type_filter)
    sample_table = collate_sbs_results(sample_table = sample_table, sbs_manipulation = sbs_manipulation, seq_type_filter = seq_type_filter)
    sample_table = collate_qc_results(sample_table = sample_table, seq_type_filter = seq_type_filter)
  }
  if(write_to_file){
    #write results from "slow option" to new cached results file
    write_tsv(sample_table, file = output_file)
  }
  #convenience columns bringing together related information
  if(join_with_full_metadata){
    if(!missing(these_samples_metadata)){
      meta_data = these_samples_metadata
    }else{
      meta_data = get_gambl_metadata(seq_type_filter = seq_type_filter)
    }

    full_table = left_join(meta_data, sample_table)

    #check for missing columns, add any missing columns and fill with NA
    col_check = c("ashm_MYC", "manta_MYC_sv", "ICGC_MYC_sv", "myc_ba", "ashm_BCL2", "manta_BCL2_sv", "ICGC_BCL2_sv", "bcl2_ba") #create a vector of the columns to check for
    missing = setdiff(col_check, names(full_table)) #return the missing columns
    full_table[missing] = NA #add missing columns to full_table and fill such columns with NAs

    full_table = full_table %>%
      mutate("MYC_SV_any" = case_when(ashm_MYC > 3 ~ "POS", manta_MYC_sv == "POS" ~ "POS", ICGC_MYC_sv == "POS" ~ "POS", myc_ba == "POS" ~ "POS", TRUE ~ "NEG"))

    full_table = full_table %>%
      mutate("BCL2_SV_any" = case_when(ashm_BCL2 > 3 ~ "POS", manta_BCL2_sv == "POS" ~ "POS", ICGC_BCL2_sv == "POS" ~ "POS", bcl2_ba == "POS" ~ "POS", TRUE ~ "NEG"))

    full_table = full_table %>%
      mutate("DoubleHitBCL2" = ifelse(BCL2_SV_any == "POS" & MYC_SV_any == "POS", "Yes", "No"))
    return(full_table)
  }
  return(sample_table)
}


#' @title Collate Derived Results.
#'
#' @description Extract derived results stored in the database (these are usually slower to derive on the fly).
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations. Default is FALSE, TRUE is not yet implemented.
#'
#' @return Data frame with one row per sample. Contains the contents of the derived_data table in the database.
#'
#' @import dplyr DBI RMariaDB
#' 
#' @noRd
#'
#' @examples
#' gambl_results_derived = collate_derived_results(samples_df)
#' 
collate_derived_results = function(sample_table,
                                   seq_type_filter = "genome",
                                   from_flatfile = FALSE){

  if(from_flatfile){
    message("not implemented YET")
  }else{
    database_name = check_config_value(config::get("database_name"))
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    derived_tbl = dplyr::tbl(con, "derived_data") %>%
      as.data.frame()
  }
  derived_tbl = derived_tbl %>%
    dplyr::select(where( ~!all(is.na(.x)))) %>%
    as.data.frame() #drop the columns that are completely empty

  print(derived_tbl)
  sample_table = dplyr::left_join(sample_table, derived_tbl)
  print(sample_table)
  return(sample_table)
}


#' @title Collate CSR Results
#'
#' @description Collate a few CSR annotations, including MiXCR.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return The sample table with additional columns.
#'
#' @import readr dplyr
#' 
#' @noRd
#'
#' @examples
#' gambl_results_derived = collate_csr_results(gambl_results_derived)
#' 
collate_csr_results = function(sample_table,
                               seq_type_filter = "genome"){

   if(seq_type_filter=="capture"){
     return(sample_table) #result doesn't exist or make sense for this seq_type
   }
   csr = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-nthomas/results/icgc_dart/mixcr_current/level_3/mixcr_genome_CSR_results.tsv"))
   sm_join = inner_join(sample_table, csr, by = c("sample_id" = "sample"))
   pt_join = inner_join(sample_table, csr, by = c("patient_id" = "sample"))
   complete_join = bind_rows(pt_join, sm_join) %>%
     bind_rows(dplyr::filter(sample_table, !patient_id %in% c(pt_join$patient_id, sm_join$patient_id))) %>%
     unique()

  return(complete_join)
}


#' @title Collate SSM Results.
#'
#' @description Compute summary statistics based on SSM calls.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param projection Specifies the projection, default is "grch37".
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations, default is TRUE.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE.
#'
#' @return The sample table with additional columns.
#'
#' @import dplyr
#' 
#' @noRd
#'
#' @examples
#' ssm_results = colalte_ssm_results(sample_table = samples,
#'                                   include_silent = TRUE)
#' 
collate_ssm_results = function(sample_table,
                               seq_type_filter = "genome",
                               projection = "grch37",
                               from_flatfile = TRUE,
                               include_silent = FALSE){

  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  seq_type = seq_type_filter
  #iterate over every sample and compute some summary stats from its MAF
  if(from_flatfile){
    base_path = check_config_value(config::get("project_base"))
    #test if we have permissions for the full gambl + icgc merge
    maf_partial_path = check_config_value(config::get("results_flatfiles")$ssm$template$merged$deblacklisted)

    maf_path = paste0(base_path, maf_partial_path)
    maf_path = glue::glue(maf_path)
    message(paste("Checking permissions on:",maf_path))
    maf_permissions = file.access(maf_path, 4)
    if(maf_permissions == -1){
      message("fail. You do not have permissions to access all the results. Use the cached results instead.")
      return(sample_table)
    }
    print(paste("loading",maf_path))
    muts = vroom::vroom(maf_path) %>%
      dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,t_alt_count,t_ref_count)
    mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
    message(paste("mutations from", mutated_samples, "samples"))
  }
  #get tally of total per sample
  muts = muts %>%
    dplyr::rename("sample_id" = "Tumor_Sample_Barcode")

  muts = mutate(muts, vaf = t_alt_count/(t_alt_count + t_ref_count))
  muts_count = dplyr::select(muts, sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("total_ssm" = "n")

  sample_table = left_join(sample_table, muts_count)
  muts_mean = muts %>%
    dplyr::select(sample_id, vaf) %>%
    group_by(sample_id) %>%
    summarize(mean_vaf = mean(vaf))

  coding_mut = dplyr::filter(muts, Variant_Classification %in% coding_class)
  coding_mut_count = coding_mut %>%
    dplyr::select(sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("coding_ssm" = "n")

  sample_table = left_join(sample_table, muts_mean)
  sample_table = left_join(sample_table, coding_mut_count)
  #check for coding SSMs in lymphoma genes
  coding_nhl = coding_mut %>%
    dplyr::filter(Hugo_Symbol %in% lymphoma_genes$Gene)

  coding_nhl_count = coding_nhl %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("driver_ssm" = "n")

  return(sample_table)
}


#' @title Collate Curated SV Results.
#'
#' @description Collate all SV calls from the genome data and summarize for main oncogenes of interest per sample.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return The sample table with additional columns.
#'
#' @import readr dplyr
#' 
#' @noRd
#'
#' @examples
#' gambl_results_derived = collate_curated_sv_results(gambl_results_derived)
#' 
collate_curated_sv_results = function(sample_table,
                                      seq_type_filter = "genome"){

  path_to_files = check_config_value(config::get("derived_and_curated"))
  project_base = check_config_value(config::get("project_base"))
  manual_files = dir(paste0(project_base, path_to_files), pattern = ".tsv")
  for(f in manual_files){
    full = paste0(project_base, path_to_files, f)
    this_data = suppressMessages(read_tsv(full, comment = "#"))
    #TO DO: fix this so it will join on biopsy_id or sample_id depending on which one is present, Done?
    sample_table = left_join(sample_table, this_data)
  }
  return(sample_table)
}


#' @title Get Sample Wildcards
#'
#' @description Get wildcards for a sample_id/seq_type combination.
#'
#' @details Return sample wildcards, useful for getting wildcard information necessary for retrieving sample-level flat-files with glue.
#'
#' @param this_sample_id The sample ID of interest.
#' @param seq_type The desired seq type, e.g genome/capture.
#'
#' @return Nothing.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' get_sample_wildcards(this_sample_id = "00-15201_tumorA",
#'                      seq_type = "genome")
#'
get_sample_wildcards = function(this_sample_id,
                                seq_type){

  sample_meta = get_gambl_metadata(seq_type_filter = seq_type) %>%
    dplyr::filter(sample_id==this_sample_id)
  this_patient_id = sample_meta$patient_id
  if(sample_meta$pairing_status=="matched"){
    normal_meta = get_gambl_metadata(seq_type_filter = seq_type,tissue_status_filter = c("normal","tumour")) %>%
      dplyr::filter(patient_id==this_patient_id) %>% dplyr::filter(tissue_status=="normal")
      normal_id = normal_meta$sample_id
      return(list(tumour_sample_id=this_sample_id,
               normal_sample_id=normal_id,
               seq_type = seq_type,
               pairing_status=sample_meta$pairing_status,
               genome_build=sample_meta$genome_build,
               unix_group=sample_meta$unix_group))
  }else{
    message("This function only works with matched samples for now")
    return()
  }
}


#' @title Assign CN to SSM.
#'
#' @description Annotate mutations with their copy number information.
#'
#' @details This function takes a sample ID with the `this_sample_id` parameter and annotates mutations with copy number information.
#' A variety of parameters are at hand for a customized workflow. For example, the user can specify if only coding mutations are of interest.
#' To do so, set `coding_only = TRUE`. It is also possible to point the function to already loaded maf/seq files, or a path to these files.
#' See parameters; `maf_file`, `maf_path`, `seq_file` and `seg_path` for more information on how to use these parameters.
#' This function can also take a vector with genes of interest (`genes`) that the returned data frame will be restricted to.
#' Is this function not what you are looking for? Try one of the following, similar, functions; [GAMBLR::get_cn_segments], [GAMBLR::get_cn_states], [GAMBLR::get_sample_cn_segments]
#'
#' @param this_sample_id Sample ID of the sample you want to annotate.
#' @param coding_only Optional. set to TRUE to restrict to only coding variants.
#' @param from_flatfile Optional. Instead of the database, load the data from a local MAF and seg file.
#' @param use_augmented_maf Boolean statement if to use augmented maf, default is FALSE.
#' @param tool_name name of tool to be used, default is "battenberg".
#' @param maf_file Path to maf file.
#' @param maf_df Optional. Use a maf dataframe instead of a path.
#' @param seg_file path to seq file.
#' @param seg_file_source Specify what copy number calling program the input seg file is from, as it handles ichorCNA differently than WisecondorX, Battenberg, etc.
#' @param assume_diploid Optional. If no local seg file is provided, instead of defaulting to a GAMBL sample, this parameter annotates every mutation as copy neutral.
#' @param genes Genes of interest.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE
#' @param this_seq_type Specified seq type for returned data.
#' @param projection specified genome projection that returned data is in reference to.
#'
#' @return A list containing a data frame (MAF-like format) with two extra columns:
#' log.ratio is the log ratio from the seg file (NA when no overlap was found)
#' as well as the segmented copy number data with the same copy number information
#' CN is the rounded absolute copy number estimate of the region based on log.ratio (NA when no overlap was found)
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr readr RMariaDB DBI ssh
#' @export
#'
#' @examples
#' cn_list = assign_cn_to_ssm(this_sample_id = "HTMCP-01-06-00422-01A-01D",
#'                            coding_only = TRUE)
#'
assign_cn_to_ssm = function(this_sample_id,
                            coding_only = FALSE,
                            from_flatfile = TRUE,
                            use_augmented_maf = TRUE,
                            tool_name = "battenberg",
                            maf_file,
                            maf_df,
                            seg_file,
                            seg_file_source = "battenberg",
                            assume_diploid = FALSE,
                            genes,
                            include_silent = FALSE,
                            this_seq_type = "genome",
                            projection = "grch37"){

  seq_type = this_seq_type

  remote_session = check_remote_configuration(auto_connect = TRUE)
  database_name = check_config_value(config::get("database_name"))
  project_base = check_config_value(config::get("project_base"))
  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  if(!missing(maf_file)){
    maf_sample = fread_maf(maf_file) %>%
      dplyr::mutate(Chromosome = gsub("chr", "", Chromosome))
  }
  else if(!missing(maf_df)){
    maf_sample = maf_df %>%
      dplyr::mutate(Chromosome = gsub("chr", "", Chromosome))
  }
  else if(from_flatfile){
    #get the genome_build and other wildcards for this sample
    wildcards = get_sample_wildcards(this_sample_id,seq_type)
    genome_build = wildcards$genome_build
    unix_group = wildcards$unix_group
    seq_type = wildcards$seq_type
    tumour_sample_id = wildcards$tumour_sample_id
    normal_sample_id = wildcards$normal_sample_id
    pairing_status = wildcards$pairing_status
    maf_sample = get_ssm_by_sample(this_sample_id = this_sample_id, this_seq_type = this_seq_type, augmented = use_augmented_maf)

  }else{
    #get all the segments for a sample and filter the small ones then assign CN value from the segment to all SSMs in that region
    con = dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    maf_table = check_config_value(config::get("results_tables")$ssm)
    maf_sample = dplyr::tbl(con, maf_table) %>%
      dplyr::filter(Tumor_Sample_Barcode == this_sample_id) %>%
      as.data.frame()
  }
  if(coding_only){
    maf_sample = dplyr::filter(maf_sample, Variant_Classification %in% coding_class)
  }

  if(!missing(genes)){
    maf_sample = dplyr::filter(maf_sample, Hugo_Symbol %in% genes)
  }

  if(!missing(seg_file)){
    seg_sample = suppressMessages(read_tsv(seg_file)) %>%
      dplyr::mutate(size = end - start) %>%
      dplyr::filter(size > 100)

    colnames(seg_sample)[c(1:4)] = c("ID", "chrom", "start", "end")
    seg_sample = seg_sample %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
      dplyr::rename(Chromosome = chrom, Start_Position = start, End_Position = end) %>%
      data.table::as.data.table()

    data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    a = data.table::as.data.table(maf_sample)
  }else if(assume_diploid == TRUE){
    if(missing(seg_file)){
      print("WARNING: A seg file was not provided! Annotating all mutation calls as copy neutral")
    }

    a = data.table::as.data.table(maf_sample)
    a_diploid = dplyr::mutate(a, CN = 2)
    return(list(maf = a_diploid))

  }else if(from_flatfile){
    message(paste("trying to find output from:", tool_name))
    project_base = check_config_value(config::get("project_base",config="default"))
    local_project_base = check_config_value(config::get("project_base"))

    results_path_template = check_config_value(config::get("results_flatfiles")$cnv$battenberg)
    results_path = paste0(project_base, results_path_template)
    local_results_path = paste0(local_project_base, results_path_template)

    ## NEED TO FIX THIS TO contain tumour/normal ID from metadata and pairing status
    battenberg_file = glue::glue(results_path)
    local_battenberg_file = glue::glue(local_results_path)


    message(paste("looking for flatfile:", battenberg_file))
    if(remote_session){
      print(local_battenberg_file)
      dirN = dirname(local_battenberg_file)

      suppressMessages(suppressWarnings(dir.create(dirN,recursive = T)))
      if(!file.exists(local_battenberg_file)){

        ssh::scp_download(ssh_session,battenberg_file,dirN)
      }
      battenberg_file = local_battenberg_file
    }

    #check for missingness
    if(!file.exists(battenberg_file)){
      print(paste("missing: ", battenberg_file))
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    seg_sample = suppressMessages(read_tsv(battenberg_file)) %>%
      as.data.table() %>%
      dplyr::mutate(size = end - start) %>%
      dplyr::filter(size > 100) %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
      dplyr::rename(Chromosome = chrom, Start_Position = start, End_Position = end)

    data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    a = data.table::as.data.table(maf_sample)
  }else{
    seg_sample = get_sample_cn_segments(this_sample_id = this_sample_id) %>%
      dplyr::mutate(size = end - start) %>%
      dplyr::filter(size > 100) %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
      dplyr::rename(Chromosome = chrom, Start_Position = start, End_Position = end) %>%
      data.table::as.data.table()

    data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    a = data.table::as.data.table(maf_sample)
    a.seg = data.table::foverlaps(a, seg_sample, type = "any")
    a$log.ratio = a.seg$log.ratio
    a$LOH = factor(a.seg$LOH_flag)
    a = dplyr::mutate(a, CN = round(2*2^log.ratio))
    seg_sample = dplyr::mutate(seg_sample, CN = round(2*2^log.ratio))
    seg_sample$LOH_flag = factor(seg_sample$LOH_flag)
  }
  if(seg_file_source == "ichorCNA"){
      message("defaulting to ichorCNA format")
      seg_sample = dplyr::rename(seg_sample, c("log.ratio" = "median", "CN" = "copy.number"))
      a.seg = data.table::foverlaps(a, seg_sample, type = "any")
      a$log.ratio = a.seg$log.ratio
      a$LOH = factor(a.seg$LOH_flag)
      a$CN = a.seg$CN

    }else{
      a.seg = data.table::foverlaps(a, seg_sample, type = "any")
      a$log.ratio = a.seg$log.ratio
      a$LOH = factor(a.seg$LOH_flag)
      a = dplyr::mutate(a, CN = round(2*2^log.ratio))
      seg_sample = dplyr::mutate(seg_sample, CN = round(2*2^log.ratio))
      seg_sample$LOH_flag = factor(seg_sample$LOH_flag)
  }
  if(!missing(seg_sample)){
    return(list(maf = a, seg = seg_sample))
  }
}


#' @title Estimate Purity.
#'
#' @description Annotate a MAF with segmented absolute copy number data and added additional columns (VAF, Ploidy and Final_purity).
#'
#' @details This function takes a sample ID with the `this_sample_id` parameter and calls [GAMBLR::assign_cn_to_ssm] to get CN information.
#' The user can also use an already loaded maf file with `maf_df`. In addition, a path to the maf/seq file of interest can also be passed to this function with
#' `in_maf` and `in_seg`. To visualize VAF and purity distributions, set the `show_plots` to TRUE (default is FALSE).
#' For more information on how to run this function with the parameters at hand, refer to the parameter descriptions and function examples.
#'
#' @param in_maf Path to a local maf file.
#' @param maf_df Optional. Instead of using the path to a maf file, use a local dataframe as the maf file.
#' @param in_seg Path to a local corresponding seg file for the same sample ID as the input maf.
#' @param this_sample_id Specify the sample_id or any other string you want embedded in the file name.
#' @param seg_file_source Specify what copy number calling program the input seg file is from, as it handles ichorCNA differently than WisecondorX, battenberg, etc.
#' @param show_plots Optional. Show two faceted plots that display the VAF and purity distributions for each copy number state in the sample. Default is FALSE.
#' @param assume_diploid Optional. If no local seg file is provided, instead of defaulting to a GAMBL sample, this parameter annotates every mutation as copy neutral. Default is FALSE.
#' @param coding_only Optional. set to TRUE to restrict to only coding variants. Default is FALSE.
#' @param genes Genes of interest.
#'
#' @return A list containing a data frame (MAF-like format) with the segmented absolute copy number data and three extra columns:
#' VAF is the variant allele frequency calculated from the t_ref_count and t_alt_count
#' Ploidy is the number of copies of an allele in the tumour cell
#' Final_purity is the finalized purity estimation per mutation after considering different copy number states and LOH events.
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' #load a maf
#' this_maf = get_ssm_by_sample(this_sample_id = "HTMCP-01-06-00422-01A-01D", 
#'                              this_seq_type = "genome")
#'
#' #estimate purity based on an already loaded maf object
#' estimate_purity(maf_df = this_maf, 
#'                 show_plots = TRUE)
#'
#' #estimate purity based sole on a smaple ID + added seg data.
#' estimate_purity(this_sample_id = "HTMCP-01-06-00422-01A-01D", 
#'                 show_plots = TRUE, coding_only = TRUE)
#'
estimate_purity = function(in_maf,
                           maf_df,
                           in_seg,
                           this_sample_id,
                           seg_file_source = "battenberg",
                           show_plots = FALSE,
                           assume_diploid = FALSE,
                           coding_only = FALSE,
                           genes){

  # Merge the CN info to the corresponding MAF file, uses GAMBLR function
  if(missing(in_maf) & missing(in_seg) & missing(maf_df)){
    CN_new = assign_cn_to_ssm(this_sample_id = this_sample_id, coding_only = coding_only, assume_diploid = assume_diploid, genes = genes, seg_file_source = seg_file_source)$maf
  }else if(!missing(in_seg)){
    CN_new = assign_cn_to_ssm(this_sample_id = this_sample_id, maf_file = in_maf, maf_df = maf_df, seg_file = in_seg, seg_file_source = seg_file_source, coding_only = coding_only, genes = genes)$maf
  }else{
    # If no seg file was provided, assume_diploid parameter is automatically set to true
    if(missing(in_seg)){
      CN_new = assign_cn_to_ssm(this_sample_id = this_sample_id, maf_file = in_maf, maf_df = maf_df, assume_diploid = TRUE, coding_only = coding_only, genes = genes)$maf
    }
  }
  # Change any homozygous deletions (CN = 0) to 1 for calculation purposes
  CN_new$CN[CN_new$CN<1] = 1
  # Select only the relevant columns, add new columns for VAF, Ploidy, and Purity
  CN_new = CN_new %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, t_ref_count, t_alt_count, CN) %>%
    dplyr::mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
    dplyr::mutate(Ploidy = ifelse(CN %in% 1:2, 1, 0)) %>%
    dplyr::mutate(Purity = "") %>%
    tidyr::drop_na(CN) %>%
    tidyr::unite("Chrom_pos", Chromosome:Start_Position, sep = ":", remove = FALSE)

  # Create an empty list
  indiv_CN = list()

  # For each copy number state, separate mutations in each CN state into separate datatables
    ## "for" loop, and this function to get different copy number states: unique(CN_new$CN)
  # Duplicate rows based on the copy number value (CN of 2 = 2 rows, etc):
    ##CN_maf[rep(seq(nrow(CN_maf)), CN_maf$CN),]
  # Fill out the ploidy column for each duplicated column:
    ##rep(seq(i),nrow(CN_max))
  # Add a column for Purity and Calculate for CNs 3 or larger:
    ##mutate(CN_max_dup, Purity = (CN*VAF)/Ploidy)

  for(i in unique(CN_new$CN)){
    CN_max = CN_new %>%
      dplyr::filter(CN == i)
    if(i > 2){
      CN_max_dup = CN_max[rep(seq(nrow(CN_max)), CN_max$CN),]
      CN_max_dup$Ploidy = rep(seq(i),nrow(CN_max))
      CN_max_dup = dplyr::mutate(CN_max_dup, Purity = (CN*VAF)/Ploidy)
      indiv_CN[[as.character(i)]] = CN_max_dup
    }else(indiv_CN[[as.character(i)]] = CN_max)
  }

  # Merge all copy state tables together into one table
  merged_CN = do.call("rbind", indiv_CN)

  # Filter for copy nnumber states 1 or 2
  # Calculate purity, and if the number is larger 1 (100%) use the VAF calue instead
  merged_CN_neut = merged_CN %>%
    dplyr::filter(merged_CN$CN < 3) %>%
    dplyr::mutate(Purity = (CN*VAF)/Ploidy) #%>%
    #dplyr::mutate(Purity = ifelse(Purity > 1, VAF, Purity))

    # Calculate a temporary purity based on the mean of these purity values
    merged_CN_neut <- merged_CN_neut %>%
      drop_na(Purity)

    mean_neut_purity = mean(merged_CN_neut$Purity)

  # For CN of 3 or larger:
    ## Filter for copy number states 1 or 2
    ## Calculate purity, and if the number is larger 1 (100%) use the VAF calue instead
    ## Group by chromosonal position,
    ## For each group, choose only the purity value that is closest to the mean_neut_purity (the temporary purity calculated from the copy neutral values)
  if (is.na(mean_neut_purity)){
    merged_CN_gain = merged_CN %>%
      dplyr::filter(merged_CN$CN > 2) %>%
      dplyr::mutate(Purity = (VAF*2)) %>%
      #dplyr::mutate(Purity = ifelse(Purity > 1, VAF, Purity)) %>%
      group_by(Chrom_pos) %>%
      slice_head() %>%
      dplyr::mutate(Assumed_CN = 2) %>%
      dplyr::mutate(Assumed_ploidy = 1)
  }else{
    merged_CN_gain = merged_CN %>%
      dplyr::filter(merged_CN$CN > 2) %>%
      dplyr::mutate(Purity = (CN*VAF)/Ploidy) %>%
      #dplyr::mutate(Purity = ifelse(Purity > 1, VAF, Purity)) %>%
      group_by(Chrom_pos) %>%
      slice(which.min(abs(Purity - mean_neut_purity)))
  }

  # Bind both datatables back together (the first contains CNs 1 and/or 2, the second contains CNs 3 and/or higher)
  CN_final = bind_rows(merged_CN_neut, merged_CN_gain)

  # Estimate the mean of all purity values from all available copy number states
  sample_purity_estimation = mean(CN_final$Purity)

  # If the final sample purity is above 1, make it equal to 1
  if(sample_purity_estimation > 1){
    sample_purity_estimation = 1
  }

  # Calculate CCF (cancer cell fraction) 2 ways:
   ## With the maximum purity estimation for the mutations in the same (largest value will be 1 for the mutation with the highest purity estimation)
   ## With the purity estimation of the sample in total. This will give values over 1 though which is problematic, as the maximum value should be 1 (100%) for clonal events
  CN_final = CN_final %>%
    dplyr::mutate(CCF_mutation = Purity/max(Purity)) %>%
    dplyr::mutate(CCF_sample = Purity/sample_purity_estimation)

  output = list()
  if(show_plots){
    # Figure 1 : VAF distribution
    # Creates facet wraps showing the VAF distribution of CN-annotated mutations for each available copy number state
    VAF_plot = CN_final %>%
      ggplot(aes(x = VAF)) +
      geom_histogram() +
      facet_wrap(~CN)

    # Figure 2: Final purity distribution
    Purity_plot = CN_final %>%
      ggplot(aes(x = Purity)) +
      geom_histogram() +
      facet_wrap(~CN)

     output[["VAF_plot"]] = VAF_plot
     output[["Purity_plot"]] = Purity_plot
  }
  output[["sample_purity_estimation"]] = sample_purity_estimation
  output[["CN_final"]] = CN_final
  return(output)
}


#' @title Refresh Full Table
#'
#' @description Refresh the contents of a database table.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::referesh_metadata_tables], not meant for out-of-package usage.
#'
#' @param table_name Name of table to refresh.
#' @param connection Database connection object.
#' @param file_path Path to the table contents to populate.
#'
#' @return A table.
#'
#' @import DBI RMariaDB readr
#' 
#' @noRd
#'
#' @examples
#' refresh_full_table(table_x, con,file_x)
#' 
refresh_full_table = function(table_name,
                              connection,
                              file_path){

  table_data = suppressMessages(read_tsv(file_path))
  dbWriteTable(con, table_name, table_data, overwrite = TRUE)
  print(paste("POPULATING table:", table_name, "USING path:", file_path))
}


#' @title Refresh Metadata Tables
#'
#' @description Refresh the contents of a metadata table.
#'
#' @details INTERNAL FUNCTION, not meant for out-of-package usage.
#'
#' @return Table.
#'
#' @import RMariaDB DBI dplyr
#' 
#' @noRd
#'
#' @examples
#' ref_meta = referesh_metadata_tables()
#' 
referesh_metadata_tables = function(){

  con = dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  all_metadata_info = sanity_check_metadata()
  tables = pull(all_metadata_info, table)
  files = pull(all_metadata_info, file)
  for(i in c(1:length(files))){
    refresh_full_table(tables[i], con, files[i])
  }
}


#' @title Sanity Check Metadata.
#'
#' @description Function that performs sanity checks on metadata.
#'
#' @details Helper function for sanity checking GAMBL metadata.
#'
#' @return A table.
#'
#' @import tibble readr dplyr
#' 
#' @noRd
#'
#' @examples
#' sane_meta_data = sanity_check_metadata()
#' 
sanity_check_metadata = function(){

  cfg = check_config_value(config::get("tables"))
  database_name = check_config_value(config::get("database_name"))
  metadata_tables = tibble(key = names(cfg), table = cfg) %>%
    unnest_auto("table")

  cfg = check_config_value(config::get("table_flatfiles"))
  metadata_files = tibble(key = names(cfg), file = cfg) %>%
    unnest_auto("file")

  all_metadata_info = left_join(metadata_tables, metadata_files)
  base_path = check_config_value(config::get("repo_base"))
  all_metadata_info = all_metadata_info %>%
    mutate(file = paste0(base_path, file))

  all_metadata_df = all_metadata_info %>%
    column_to_rownames(var = "key")
  #all samples with different seq_type and protocol must have a unique sample_id
  sample_df = suppressMessages(read_tsv(all_metadata_df["samples", "file"]))
  tumour_samples = sample_df %>%
    dplyr::select(patient_id, sample_id, biopsy_id, seq_type, protocol) %>%
    dplyr::filter(!is.na(biopsy_id))

  n_samp_bio = tumour_samples %>%
    count() %>%
    pull(n)

  #check for any multiplicity of sample_id
  n_samp = tumour_samples %>%
    dplyr::select(-biopsy_id) %>%
    unique() %>%
    count() %>%
    pull(n)

  #should be the same number as above
  if(!n_samp == n_samp_bio){
    print("ERROR! some biopsies appear to have the same sample_id/protocol combination")
  }
  return(all_metadata_info)
}


#' @title Collate Ancestry.
#'
#' @description Gather ancestry information and expand the incoming sample table (or metadata).
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param somalier_output Somalier ancestery.tsv
#'
#' @return A table.
#'
#' @import stringr readr dplyr
#' 
#' @noRd
#'
#' @examples
#' table = collate_ancestry(sample_table = "my_sample_table.txt")
#' 
collate_ancestry = function(sample_table,
                            seq_type_filter="genome",
                            somalier_output){
  if(seq_type_filter=="capture"){
    message("skipping ancestry for this seq_type")
    return(sample_table)
  }
  if(missing(somalier_output)){
    somalier_output = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/gambl/somalier_current/02-ancestry/2020_08_07.somalier-ancestry.tsv"
  }
  somalier_all = suppressMessages(read_tsv(somalier_output))
  somalier_all = mutate(somalier_all, sample_id = str_remove(`#sample_id`, pattern = ":.+")) %>%
    dplyr::select(-`#sample_id`, -given_ancestry)
  somalier_all = dplyr::select(somalier_all, sample_id, predicted_ancestry, PC1, PC2, PC3, PC4, PC5)
  sample_table = left_join(sample_table, somalier_all)
  return(sample_table)
}


#' @title Collate Extra Metadata.
#'
#' @description Gather additional metadata information and expand the incoming sample table (or metadata).
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param file_path Path to extra metadata.
#'
#' @return A table.
#'
#' @import readr dplyr
#' 
#' @noRd
#'
#' @examples
#' table = collate_extra_metadata(sample_table = "my_sample_table.txt")
#' 
collate_extra_metadata = function(sample_table,
                                  file_path){

  file_path = "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/experiments/2021-04-21-Trios-MiXCR/trios_relapse_timing.tsv"
  extra_df = suppressMessages(read_tsv(file_path))
  sample_table = left_join(sample_table, extra_df, by = c("sample_id" = "biopsy_id"))
}


#' @title Collate SBS Results.
#'
#' @description Bring in the results from mutational signature analysis.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param file_path Optional path to SBS file.
#' @param scale_vals Parameter not used?
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#'
#' @return A data frame with new columns added.
#'
#' @import dplyr tibble
#' 
#' @noRd
#'
#' @examples
#' collated = collate_sbs_results(sample_table = sample_table,
#'                                sbs_manipulation = sbs_manipulation)
#
collate_sbs_results = function(sample_table,
                               seq_type_filter = "genome",
                               file_path,
                               scale_vals = FALSE,
                               sbs_manipulation = ""){
  if(seq_type_filter!="genome"){
    message("skipping sbs for seq_type")
    return(sample_table)
  }
  if(missing(file_path)){
    base = check_config_value(config::get("project_base"))

    file_path = paste0(base,"icgc_dart/sigprofiler-1.0/02-extract/genome--hg38/BL_HGBL_DLBCL_FL_COMFL_CLL_MCL_B-ALL_PBL_DLBCL-BL-like_UNSPECIFIED_SCBC_MM_all/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt")
  }
  message(paste("loading",file_path))
  signatures = read.csv(file_path, sep = "\t", header = 1, row.names = 1)
  rs = rowSums(signatures)
  cn = colnames(signatures)
  new_sig = signatures
  if(sbs_manipulation == "scale"){
    #for(col in cn){
    #  scaled_vals = signatures[,col] / rs
    #  new_sig[,col]=scaled_vals
    #}
    #sbs1 = signatures[,"SBS1"] / rs
    #sbs5 = signatures[,"SBS5"] / rs
    #sbs9 = signatures[,"SBS9"] / rs
    #sbs8 = signatures[,"SBS8"] / rs
    #sbs12 = signatures[,"SBS12"] / rs
    #sbs17b = signatures[,"SBS17b"]/rs
    #sbs18 = signatures[,"SBS18"]/rs
    #sbs84 = signatures[,"SBS84"]/rs
    #sbs85 = signatures[,"SBS85"]/rs
    #sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs5=sbs5,sbs9=sbs9,sbs8=sbs8,sbs12=sbs12,sbs17b=sbs17b,sbs18=sbs18,sbs84=sbs84,sbs85=sbs85)
    sbs = apply(signatures, 2, function(x){x/rs}) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
  }
  else if(sbs_manipulation == "log"){
    #sbs1 = log(signatures[,"SBS1"]+1)
    #sbs5 = log(signatures[,"SBS5"]+1)
    #sbs9 = log(signatures[,"SBS9"]+1)
    #sbs8 = log(signatures[,"SBS8"]+1)
    #sbs12 = log(signatures[,"SBS12"]+1)
    #sbs17b = log(signatures[,"SBS17b"]+1)
    #sbs18 = log(signatures[,"SBS18"]+1)
    #sbs84 = log(signatures[,"SBS84"]+1)
    #sbs85 = log(signatures[,"SBS85"]+1)
    sbs = apply(signatures, 2, function(x){log(x + 1)}) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
    #sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs5=sbs5,sbs9=sbs9,sbs8=sbs8,sbs12=sbs12,sbs17b=sbs17b,sbs18=sbs18,sbs84=sbs84,sbs85=sbs85)
  }else{
    #sbs1 = signatures[,"SBS1"]
    #sbs5 = signatures[,"SBS5"]
    #sbs9 = signatures[,"SBS9"]
    #sbs8 = signatures[,"SBS8"]
    #sbs12 = signatures[,"SBS12"]
    #sbs17b = signatures[,"SBS17b"]
    #sbs18 = signatures[,"SBS18"]
    #sbs84 = signatures[,"SBS84"]
    #sbs85 = signatures[,"SBS85"]
    #sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs5=sbs5,sbs9=sbs9,sbs8=sbs8,sbs12=sbs12,sbs17b=sbs17b,sbs18=sbs18,sbs84=sbs84,sbs85=sbs85)
    sbs = signatures %>%
      rownames_to_column("sample_id")
  }
  sample_table = left_join(sample_table, sbs)
  return(sample_table)
}


#' @title Collate NFKBIZ Results.
#'
#' @description Determine which cases have NFKBIZ UTR mutations.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return Samples table.
#'
#' @import dplyr
#' 
#' @noRd
#'
#' @examples
#' sample_table = collate_nfkbiz_results(sample_table = sample_table)
#' 
collate_nfkbiz_results = function(sample_table,
                                  seq_type_filter = "genome"){

  #TO DO: Update to work with hg38 projection
  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter=seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  this_region = "chr3:101578214-101578365"
  nfkbiz_ssm = get_ssm_by_region(region = this_region,seq_type = seq_type_filter) %>%
    pull(Tumor_Sample_Barcode) %>%
    unique
  if(seq_type_filter=="genome"){
    nfkbiz_sv = get_manta_sv(region = this_region) %>%
      pull(tumour_sample_id) %>%
      unique
    nfkbiz = unique(c(nfkbiz_ssm, nfkbiz_sv))
  }
  else{
    nfkbiz = unique(nfkbiz_ssm)
  }


  sample_table$NFKBIZ_UTR = "NEG"
  sample_table[sample_table$sample_id %in% nfkbiz, "NFKBIZ_UTR"] = "POS"
  return(sample_table)
}


#' @title Collate ASHM Results.
#'
#' @description Determine the hypermutation status of a few genes.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return Samples table.
#'
#' @import dplyr tidyr tibble
#' 
#' @noRd
#'
#' @examples
#' sample_table = collate_ashm_results(sample_table = sample_table)
#' 
collate_ashm_results = function(sample_table,
                                seq_type_filter = "genome"){

  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter = seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  #just annotate BCL2, MYC and CCND1 hypermutation
  regions_df = data.frame(name = c("CCND1","BCL2","MYC"),
  region = c("chr11:69455000-69459900", "chr18:60983000-60989000", "chr8:128747615-128751834"))
  region_mafs = lapply(regions_df$region, function(x){get_ssm_by_region(region = x, streamlined = FALSE,seq_type=seq_type_filter)})
  tibbled_data = tibble(region_mafs, region_name = regions_df$name)
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)

  unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
    dplyr::select(start, sample_id, region_name)

  tallied = unlisted_df %>%
    group_by(sample_id, region_name) %>%
    tally() %>%
    pivot_wider(values_from = n, names_from = region_name, values_fill = 0, names_prefix = "ashm_")

  sample_table = left_join(sample_table, tallied, by = "sample_id")
}


#' @title Collate SV Results.
#'
#' @description Determine and summarize which cases have specific oncogene SVs.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param tool Name of tool (optional, default is manta).
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param oncogenes Which oncogenes to collate SVs from.
#'
#' @return Data frame with additional columns ({tool}_{oncogene} and {tool}_{oncogene}_{partner}).
#'
#' @import dplyr
#' 
#' @noRd
#'
#' @examples
#' results = collate_samples_sv_results(sample_table = samples,
#'                                      tool = "manta",
#'                                      oncogenes = c("MYC", "BCL2"))
#' 
collate_sv_results = function(sample_table,
                              tool = "manta",
                              seq_type_filter = "genome",
                              oncogenes = c("MYC", "BCL2", "BCL6", "CCND1", "IRF4")){
  if(seq_type_filter!="genome"){
    message("skipping sv for this seq_type")
    return(sample_table)
  }
  if(tool == "manta"){
    all_svs = get_manta_sv()
  }
  annotated_svs = annotate_sv(all_svs) %>%
  dplyr::filter(!is.na(partner))

  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  multiout = function(df,
                       annotated,
                       tool,
                       oncogene_name){

    some_fusions = dplyr::filter(annotated, gene == all_of(oncogene_name)) %>%
      group_by(tumour_sample_id) %>%
      arrange(partner) %>%
      dplyr::filter(row_number() == 1)

    df = mutate(df, "{tool}_{oncogene_name}_sv" := case_when(sample_id %in% some_fusions$tumour_sample_id ~ "POS", TRUE ~ "NEG"))
    some_fusions = some_fusions %>%
      dplyr::select(tumour_sample_id, partner) %>%
      mutate("{tool}_{oncogene_name}_partner" := partner) %>%
      dplyr::select(-partner)

    df = left_join(df, some_fusions, by = c("sample_id" = "tumour_sample_id"))
    return(df)
  }
  out_table = sample_table
  for(oncogene in oncogenes){
    out_table = multiout(out_table, annotated_svs, "manta", oncogene)
  }
  return(out_table)
}


#' @title Get GAMBL Colours.
#'
#' @description Get GAMBL colour schemes for annotating figures.
#'
#' @details This function was designed to retrieve specified GAMBL colour palettes.
#' By default, this function returns all the colours currently available.
#' The user can easily specify what classification to return colors for with the `classification` parameter.
#' It is also possible to return any given colour in different formats.
#' To do so, refer to the Boolean arguments; `as_list` and `as_dataframe`.
#' For more information regarding the available colours, refer to the utilities vignette.
#'
#' @param classification Optionally request only colours for pathology, lymphgen, mutation or copy_number.
#' @param alpha Alpha of plotted colours.
#' @param as_list Boolean parameter controlling the format of the return. Default is FALSE.
#' @param as_dataframe Boolean parameter controlling the format of the return. Default is FALSE.
#' @param return_available Set to TRUE for returning all available colours. Default is FALSE.
#' @param verbose Default is FALSE
#'
#' @return A named vector of colour codes for lymphgen classes and pathology.
#'
#' @import dplyr ggsci stringr tidyr
#' @export
#'
#' @examples
#' lymphgen_cols = get_gambl_colours("lymphgen")
#' 
#' \dontrun{
#' #be sure to install ggsci from https://github.com/morinlab/ggsci
#' #install_github("morinlab/ggsci")
#' }
#' 
get_gambl_colours = function(classification = "all",
                             alpha = 1,
                             as_list = FALSE,
                             as_dataframe = FALSE,
                             return_available = FALSE,
                             verbose = FALSE){

  all_colours = list()
  everything = c()
  blood_cols = ggsci::get_ash("blood")

  all_colours[["seq_type"]] = c("mrna" = "#E41A1C",
                                "genome" = "#377EB8",
                                "capture" = "#4DAF4A")

  all_colours[["type"]] = c("gain" = "blue",
                            "loss" = "red")

  all_colours[["hmrn"]] = c("BCL2-MYC" = "#52000F",
                            "BCL2" = "#721F0F",
                            "SOCS1/SGK1" = "#D66B1F",
                            "TET2/SGK1" = "#C41230",
                            "MYD88" = "#3B5FAC",
                            "NOTCH2" = "#7F3293",
                            "NOTCH1" = "#55B55E",
                            "Other" = "#ACADAF")

  all_colours[["EBV"]] =  c("EBV-positive" = "#7F055F",
                            "EBV-negative" = "#E5A4CB",
                            "POS" = "#7F055F",
                            "NEG" = "#E5A4CB")

  all_colours[["BL"]] = c("Q53-BL" = "#A6CEE3",
                          "M53-BL" = "#A6CEE3", #added because genetic subgroup still refers to it this way
                          "DLBCL-A" = "#721F0F",
                          "IC-BL" = "#45425A",
                          "DGG-BL" = "#E90C8B",
                          "DLBCL-B" = "#FB9A99",
                          "DLBCL-C" = "#C41230")

  all_colours[["FL"]] = c(dFL = "#99C1B9", cFL = "#D16666", DLBCL = "#479450")

  all_colours[["lymphgenerator"]] = c("MP3"="#5B8565",
                                      "EGB" = "#98622A",
                                      "ETB"="#813F3D",
                                      "aSCI"="#D66B1F",
                                      "aSEL"="#6A0D18",
                                      "MCaP"="#5F8CFF",
                                      "BNZ"="#8870B6",
                                      "EZB"="#721F0F",
                                      "ST2"="#C41230",
                                      "UNCLASS"="#05631E"
                                      )

  all_colours[["chapuy_classifier"]] = c(
    C0 = "#bebebe",
    C1 = "#803D99",
    C2 ="#00A2D2",
    C3 = "#F39123",
    C4 = "#50BFAD",
    C5 = "#DE292A"
  )

  all_colours[["lacy_classifier"]] = all_colours[["hmrn"]]

  all_colours[["lymphgen"]] = c("EZB-MYC" = "#52000F",
                                "EZB" = "#721F0F",
                                "EZB-COMP" = "#C7371A",
                                "ST2" = "#C41230",
                                "ST2-COMP" = "#EC3251",
                                "MCD" = "#3B5FAC",
                                "MCD-COMP" = "#6787CB",
                                "BN2" =  "#7F3293",
                                "BN2-COMP" = "#A949C1",
                                "N1" = "#55B55E",
                                "N1-COMP" = "#7FC787",
                                "A53" = "#5b6d8a",
                                "Other" = "#ACADAF",
                                "COMPOSITE" = "#ACADAF")

  #all_colours[["coding_class"]] = c("Frame_Shift_Del","Frame_Shift_Ins",
  #                 "In_Frame_Del","In_Frame_Ins",
  #                 "Missense_Mutation","Nonsense_Mutation",
  #                 "Nonstop_Mutation","Splice_Region","Splice_Site",
  #                 "Targeted_Region","Translation_Start_Site")
  all_colours[["mutation"]]=
    c(
        "Nonsense_Mutation"=unname(blood_cols["Red"]),
        "Missense_Mutation"=unname(blood_cols["Green"]),
        "Multi_Hit"=unname(blood_cols["Steel Blue"]),
        "Frame_Shift_Ins" = unname(blood_cols["Magenta"]),
        "Frame_Shift_Del" = unname(blood_cols["Magenta"]),
        "In_Frame_Ins" = unname(blood_cols["Brown"]),
        "In_Frame_Del" = unname(blood_cols["Brown"]),
        "Nonstop_Mutation" = unname(blood_cols["Light Blue"]),
        "Translation_Start_Site" = unname(blood_cols["Lavendar"]),
        "Splice_Site" = unname(blood_cols["Orange"]),
        "Splice_Region" = unname(blood_cols["Orange"]),
        "3'UTR" = unname(blood_cols["Yellow"]),
        "Silent" = "#A020F0")

  all_colours[["rainfall"]] =
    c(
      "C>A" = "#2196F3FF",
      "C>G" = "#3F51B5FF",
      "C>T" = "#F44336FF",
      "InDel" = "purple",
      "T>A" = "#4CAF50FF",
      "T>C" = "#FFC107FF",
      "T>G" = "#FF9800FF"
    )

  all_colours[["pos_neg"]]=c(
    "POS"="#c41230",
    "NEG"="#E88873",
    "PARTIAL"="#E88873",
    "yes"="#c41230",
    "no"="#E88873",
    "YES"="#c41230",
    "NO"="#E88873",
    "FAIL"="#bdbdc1",
    "positive"="#c41230",
    "negative"="#E88873",
    "fail"="#bdbdc1")

  all_colours[["copy_number"]]=c(
    "nLOH"="#E026D7",
    "14"="#380015",
    "15"="#380015",
    "13"="#380015",
    "12"="#380015",
    "11"="#380015",
    "10"="#380015",
    "9"="#380015",
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
  all_colours[["blood"]] = c(
      "Red" = "#c41230", "Blue"="#115284","Green" = "#39b54b",
      "Purple" = "#5c266c", "Orange"="#fe9003","Green" = "#046852",
      "Lavendar" = "#8781bd", "Steel Blue"= "#455564",
      "Light Blue" = "#2cace3", "Magenta" = "#e90c8b", "Mustard" = "#b76d29",
      "LimeGreen" = "#a4bb87", "Brown" = "#5f3a17", "Gray" = "#bdbdc1",
      "Yellow" = "#f9bd1f"
  )
  all_colours[["sex"]]=c(
    "M"="#118AB2",
    "Male"="#118AB2",
    "male"="#118AB2",
    "F"="#EF476F",
    "Female"="#EF476F",
    "female"="#EF476F")
  all_colours[["clinical"]]=ggsci::get_ash("clinical")
  all_colours[["pathology"]] = c(
      "B-ALL"="#C1C64B",
      "CLL"="#889BE5",
      "MCL"="#F37A20",
      "BL"="#926CAD",
      "mBL"="#34C7F4",
      "tFL"="#FF8595",
      "DLBCL-BL-like"="#34C7F4",
      "pre-HT"="#754F5B",
      "PMBL"= "#227C9D",
      "PMBCL"="#227C9D",
      "FL"="#EA8368",
      "no-HT"="#EA8368",
      "COMFL"="#8BBC98",
      "COM"="#8BBC98",
      "post-HT"="#479450",
      "DLBCL"="#479450",
      "denovo-DLBCL"="#479450",
      "HGBL-NOS"="#294936",
      "HGBL"="#294936",
      "HGBL-DH/TH"="#7A1616",
      "PBL" = "#E058C0",
      "Plasmablastic" = "#E058C0",
      "CNS" = "#E2EF60",
      "THRLBCL" = "#A5F2B3",
      "MM"="#CC9A42",
      "SCBC"="#8c9c90",
      "UNSPECIFIED"="#cfba7c",
      "OTHER"="#cfba7c",
      "MZL"="#065A7F",
      "SMZL"="#065A7F",
      "Prolymphocytic" = "#7842f5"
  )
  all_colours[["coo"]] = c(
    "ABC" = "#05ACEF",
    "UNCLASS" = "#05631E",
    "Unclass" = "#05631E",
    "U" = "#05631E",
    "UNC" = "#05631E",
    "GCB"= "#F58F20",
    "DHITsig-"= "#F58F20",
    "DHITsigNeg"= "#F58F20",
    "DHITsig-IND" = "#003049",
    "DHITsig+" = "#D62828",
    "DHITsigPos" = "#D62828",
    "NA" = "#ACADAF"
  )
  all_colours[["cohort"]] = c("Chapuy"="#8B0000","Chapuy, 2018"="#8B0000",
                  "Arthur"= "#8845A8","Arthur, 2018"= "#8845A8",
                  "Schmitz"= "#2C72B2","Schmitz, 2018"= "#2C72B2",
                  "Reddy" = "#E561C3","Reddy, 2017" = "#E561C3",
                  "Morin"= "#8DB753", "Morin, 2013"= "#8DB753",
                  "Kridel"= "#4686B7", "Kridel, 2016"= "#4686B7",
                  "ICGC"="#E09C3B","ICGC, 2018"="#E09C3B",
                  "Grande"="#e90c8b", "Grande, 2019"="#e90c8b")

  all_colours[["indels"]] = c("DEL" = "#53B1FC", "INS" = "#FC9C6D")
  all_colours[["svs"]] = c("DEL" = "#53B1FC", "DUP" = "#FC9C6D")
  all_colours[["genetic_subgroup"]] = c(all_colours[["lymphgen"]],all_colours[["BL"]],all_colours[["FL"]])
  #print(all_colours)
  if(alpha <1){
    for(colslot in names(all_colours)){
      raw_cols = all_colours[[colslot]]
      raw_cols_rgb = col2rgb(raw_cols)
      alpha_cols = rgb(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ], alpha = alpha * 255L, names = names(raw_cols), maxColorValue = 255L)
      names(alpha_cols) = names(raw_cols)
      all_colours[[colslot]] = alpha_cols
    }
  }
  for(this_group in names(all_colours)){
    everything = c(everything, all_colours[[this_group]])
  }
  #return matching value from lowercase version of the argument if it exists
  lc_class = stringr::str_to_lower(classification)
  if(return_available){
    return(names(all_colours))
  }
  if(classification %in% names(all_colours)){
    if(as_dataframe){
      some_col=all_colours[[classification]]
      df_ugly = data.frame(name=names(some_col),colour=unname(some_col))
      df_tidy = mutate(df_ugly,group=classification)
      return(df_tidy)
    }
    return(all_colours[[classification]])
  }else if(lc_class %in% names(all_colours)){
    return(all_colours[[lc_class]])
  }else if(as_list){
    return(all_colours)
  }else if(as_dataframe){
    df_ugly = data.frame(name = names(unlist(all_colours, use.names = T)), colour = unlist(all_colours, use.names = T))
    df_tidy = separate(df_ugly,name,into=c("group","name"),sep="\\.")
    return(df_tidy)
  }else{
    return(everything)
  }
}


#' @title Get BAMs.
#'
#' @description Get full paths for bam files for a sample or patient.
#'
#' @details Returns a list with BAM paths for tumour, normal and mrna data.
#' This function expects a sample ID (`this_sample_id`) or a patient ID (`this_patient_id`).
#'
#' @param this_sample_id Sample ID of interest.
#' @param this_patient_id patient ID of interest.
#'
#' @return A list that contains the genome_build and an igv-friendly build (igv_build), a list of bam file paths for tumour, normal and mrna data.
#'
#' @import dplyr
#' @export
#'
#' @examples
#'
#' #example 1, using a sample ID
#' bam_details = get_bams(this_sample_id = "HTMCP-01-06-00422-01A-01D")
#'
#' #example 2, using a patient ID
#' bam_details = get_bams(this_patient_id = "HTMCP-01-06-00422")
#'
get_bams = function(this_sample_id,
                    this_patient_id){

  meta = get_gambl_metadata(tissue_status_filter = c("tumour", "normal"), seq_type_filter = "genome")
  meta_mrna = get_gambl_metadata(seq_type_filter = "mrna")
  #get all samples for this patient
  if(missing(this_patient_id)){
    this_patient_id = meta %>%
      dplyr::filter(sample_id == this_sample_id) %>%
      dplyr::pull(patient_id)
  }
  meta_patient = meta %>%
    dplyr::filter(patient_id == this_patient_id)

  meta_mrna_patient = meta_mrna %>%
    dplyr::filter(patient_id == this_patient_id)

  build = dplyr::pull(meta_patient, genome_build) %>%
    head(1)
  if(build == "hs37d5"){
    igv_build = "hg19"
  }else{
    igv_build = build
  }
  tumour_genome_bams = dplyr::filter(meta_patient, seq_type == "genome" & tissue_status == "tumour") %>%
    dplyr::pull(data_path)

  bam_details = list(igv_build = igv_build, genome_build = build, tumour_bams = tumour_genome_bams)
  normal_genome_bams = dplyr::filter(meta_patient, seq_type == "genome" & tissue_status == "normal") %>%
    dplyr::pull(data_path)

  unix_group = dplyr::filter(meta_patient, seq_type == "genome" & tissue_status == "tumour") %>%
    dplyr::pull(unix_group) %>%
    unique()

  bam_details$pairing_status = get_gambl_metadata(seq_type_filter = "genome") %>%
    dplyr::filter(tissue_status == "tumour", patient_id == this_patient_id) %>%
    dplyr::pull(pairing_status) %>%
    unique()

  bam_details$unix_group = unix_group
  if(length(normal_genome_bams)){
    bam_details$normal_genome_bams = normal_genome_bams
  }
  if(length(normal_genome_bams)){
    bam_details$normal_genome_bams = normal_genome_bams
  }
  rnaseq_bams = dplyr::filter(meta_mrna_patient, seq_type == "mrna") %>%
    dplyr::pull(data_path)
  if(length(rnaseq_bams)){
    bam_details$rnaseq_bams = rnaseq_bams
  }
  return(bam_details)
}


#' @title Make IGV Snapshot
#'
#' @description Load bams and generate an IGV screenshot for one or more regions.
#'
#' @details Specify the path to one or more bam files as a character vector to the `bams` parameter.
#' The user can also specify regions of interest with either the `region` parameter (chr:start-end),
#' or the user can directly supply the chromosome, start and end coordinates with the `chrom`, `start`, and `end` parameters.
#' For more information and examples, refer to the function examples and parameter descriptions.
#'
#' @param bams Character vector containing the full path to one or more bam files.
#' @param genome_build String specifying the genome build for the bam files provided.
#' @param region Optionally specify the region as a single string (e.g. "chr1:1234-1235").
#' @param padding Optionally specify a positive value to broaden the region around the specified position. Default is 200.
#' @param chrom Optionally specify the region by specifying the chromosome, start and end (see below).
#' @param start Optionally specify the region by specifying the start.
#' @param end Optionally specify the region by specifying the end.
#' @param this_sample_id Specify the sample_id or any other string you want embedded in the file name.
#' @param out_path Specify the output directory where the snapshot will be written.
#' @param igv_port Specify the port IGV is listening on.
#'
#' @return Path to file (.png).
#'
#' @import SRAdb
#' @export
#'
#' @examples
#' \dontrun{
#' #IMPORTANT: you must be running IGV on the host that is running R and you need to have it
#' #listening on a port. The simplest scenario is to run this command on a terminal (if using a Mac),
#' #assuming you are using R on gphost10 and you have a ssh config that routes gp10 to that host
#' 
#' ssh -X gp10
#' 
#' #then launch IGV (e.e. from a conda installation):
#' #conda activate igv; igv &
#' this_sv = annotated_sv %>% 
#'  filter(gene=="ETV6")
#' 
#' tumour_bam = get_bams(this_sample_id = this_sv$tumour_sample_id)
#' 
#' make_igv_snapshot(chrom = this_sv$chrom2,
#'                   start = this_sv$start2,
#'                   end = this_sv$end2,
#'                   this_sample_id = this_sv$tumour_sample_id,
#'                   out_path = "~/IGV_snapshots/")
#' }
#'
make_igv_snapshot = function(bams,
                             genome_build,
                             region,
                             padding = 200,
                             chrom,
                             start,
                             end,
                             this_sample_id,
                             out_path = "/tmp/",
                             igv_port = 60506){

  sock = IGVsocket(port = igv_port)
  IGVclear(sock)
  if(missing(region)){
    region = paste0(chrom, ":", start-padding, "-", end + padding)
  }
  IGVgenome(sock, genome = genome_build)
  IGVgoto(sock, region)
  for(bam_file in bams){
    IGVload(sock, bam_file)
  }
  filename = paste(this_sample_id, region, "snapshot.png", sep = "_")
  IGVsnapshot(sock, fname = filename, dirname = out_path)
  return(paste0(out_path, filename))
}


#' @title Fishers Exact Test (CNV).
#'
#' @description Using GISTIC2.0 outputs, perform Fisher's exact test to compare CNV frequencies between 2 groups.
#'
#' @details This function was developed to compare (Fisher's exact test) CNV frequencies between two groups.
#' To do so, set the path to the GISTIC2.0 all_lesions file with `gistic_lesions`, together with a metadata table with the sample IDs of interest (`metadata`).
#' The last remaining required parameter is `comparison`, this parameter takes the name of the column annotating the groups of interest, e.g pathology, cohort, etc.
#' For more information on how to run this function, refer to the function examples and parameter descriptions.
#'
#' @param gistic_lesions Path to the GISTIC2.0 all_lesions output file.
#' @param metadata Data frame containing sample ids and column with annotated data for the 2 groups of interest. All other columns are ignored. Currently, function exists if asked to compare more than 2 groups.
#' @param comparison Specify column annotating groups of interest.
#' @param fdr.method FDR method to adjust p values. Uses p.adjust function, and therefore accepts its method for FDR ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). By default, this function uses "fdr".
#' @param fdr.cutoff Specify FDR significance cut-off. By default, this function uses 0.1.
#' @param text_size Size of the text on the forest plot of differentially enriched CNV. Default text-size is 7.
#' @param blacklisted_regions Optionally, specify any descriptors (value from column `Descriptor` of GISTIC2.0 all_lesions output file) to filter out before any comparisons are done. It is possible to specify a list of multiple descriptors, for example, c("3p12.3", "12p13.2"). Default is NULL.
#'
#' @return list
#'
#' @import dplyr metaviz readr tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' # basic usage
#' FtestCNV(gistic_lesions = "path_to_GISTIC2.0_output/all_lesions.conf_{confidence_level}.txt",
#'          metadata = derived_data,
#'          comparison = "pathology")
#' 
#' # advanced usage
#' FtestCNV(gistic_lesions = "path_to_GISTIC2.0_output/all_lesions.conf_{confidence_level}.txt",
#'          metadata = derived_data,
#'          comparison = "pathology",
#'          fdr.method = "bonferroni",
#'          fdr.cutoff = 0.05,
#'          blacklisted_regions = c("3p12.3", "12p13.2"))
#' }
#'
FtestCNV = function(gistic_lesions,
                     metadata,
                     comparison,
                     fdr.method = "fdr",
                     fdr.cutoff = 0.1,
                     text_size = 7,
                     blacklisted_regions = NULL){

  # get groups of interest for comparison
  GROUPS.TO.COMPARE = unique(metadata[,comparison]) %>% unlist

  # quick check that only 2 groups are provided forr comparison
  if(length(GROUPS.TO.COMPARE) > 2){
    message("The current implementation of function only accepts 2 groups for comparison. You provided 3 groups.")
    message("Please modify metadata accordingly to compare only 2 groups.")
    message("Groups you provided are: ", paste(c(GROUPS.TO.COMPARE), collapse = ","))
    return(NULL)
  }
  # read lesions from gistic utput to collect event/case
  lesions = suppressMessages(read_tsv(gistic_lesions, col_names = TRUE)) %>%
    dplyr::filter(!grepl("CN", `Unique Name`)) %>%
    dplyr::select(-tail(names(.), 1), -`Residual q values after removing segments shared with higher peaks`, -`Broad or Focal`, -`Amplitude Threshold`, -`q values`) %>%
    dplyr::filter (! Descriptor %in% blacklisted_regions)

  # prepare metadata
  # save names of first few columns of gistic lesions file for future reference
  columns = colnames(lesions)[1:5]
  # subset lesions file to only lesions metadata and samples of interest
  columns = c(columns, pull(metadata, Tumor_Sample_Barcode))
  lesions = subset(lesions, select = intersect(columns, colnames(lesions)))

  # subset ids for each category for future comparisons
  GROUP1 = dplyr::pull(metadata %>%
    dplyr::filter(base::get(comparison) == GROUPS.TO.COMPARE[1]), Tumor_Sample_Barcode)

  GROUP2 = dplyr::pull(metadata %>%
    dplyr::filter(base::get(comparison) == GROUPS.TO.COMPARE[2]), Tumor_Sample_Barcode)

  # subset lesions file for each group and count number of CNV
  GROUP1.CNV = lesions[,colnames(lesions) %in% GROUP1] %>%
    dplyr::mutate(count = rowSums(.!=0))

  GROUP2.CNV = lesions[,colnames(lesions) %in% GROUP2] %>%
    dplyr::mutate(count = rowSums(.!=0))

  ########## PART 1
  # compare significance of differences between 2 groups

  # first, get the counts of mutated samples for each group
  GROUP1_vs_GROUP2 = as.data.frame(cbind(lesions[,1:5], GROUP1.CNV$count, GROUP2.CNV$count))
  colnames(GROUP1_vs_GROUP2)[6:7] = c("Mutated_GROUP1", "Mutated_GROUP2")

  # second, get the counts of unmutated samples for each group
  GROUP1_vs_GROUP2 = GROUP1_vs_GROUP2 %>%
    dplyr::mutate(Not_Mutated_GROUP1 = length(GROUP1) - Mutated_GROUP1, .after = Mutated_GROUP1) %>%
    dplyr::mutate(Not_Mutated_GROUP2 = length(GROUP2) - Mutated_GROUP2, .after = Mutated_GROUP2) %>%
    dplyr::mutate(Region = ifelse(grepl("Amplification", `Unique Name`), paste0(Descriptor, "_Amp"), paste0(Descriptor, "_Del")), .after = Descriptor)

  # row-by-row, calculate 2-tailed Fisher's exact test for each CNV
  GROUP1_vs_GROUP2 = GROUP1_vs_GROUP2 %>%
    dplyr::filter(Mutated_GROUP1 / (Mutated_GROUP1 + Not_Mutated_GROUP1) > 0.05 | Mutated_GROUP2 / (Mutated_GROUP2 + Not_Mutated_GROUP2) > 0.05) %>%
    dplyr::mutate(pVal = apply(., 1, function(x) {tbl = matrix(as.numeric(x[7:10]), ncol = 2, byrow = T); fisher.test(tbl, alternative = "two.sided")$p.value})) %>%
    dplyr::mutate(OddsRatio = apply(., 1, function(x) {tbl = matrix(as.numeric(x[7:10]), ncol = 2, byrow = T); log(fisher.test(tbl, alternative = "two.sided")$estimate)})) %>%
    dplyr::mutate(LowConfInt = apply(., 1, function(x) {tbl = matrix(as.numeric(x[7:10]), ncol = 2, byrow = T); log(fisher.test(tbl, alternative = "two.sided")$conf.int[1])})) %>%
    dplyr::mutate(HighConfInt = apply(., 1, function(x) {tbl = matrix(as.numeric(x[7:10]), ncol = 2, byrow = T); log(fisher.test(tbl, alternative = "two.sided")$conf.int[2])}))

  # adjust FDR for multiple comparisons
  GROUP1_vs_GROUP2$FDR = p.adjust(GROUP1_vs_GROUP2$pVal, method = fdr.method)

  # filter for only those CNV that pass FDR cutoff
  GROUP1_vs_GROUP2.PASSED = GROUP1_vs_GROUP2[GROUP1_vs_GROUP2$FDR<fdr.cutoff,]
  GROUP1_vs_GROUP2.PASSED = GROUP1_vs_GROUP2.PASSED %>%
    dplyr::distinct(Region, .keep_all = TRUE)

  # rename columns to match with comparison groups and save resulting df as part of the output
  DISTINCT = GROUP1_vs_GROUP2.PASSED %>%
    `colnames<-`(gsub("GROUP1", GROUPS.TO.COMPARE[1], colnames(.))) %>%
    `colnames<-`(gsub("GROUP2", GROUPS.TO.COMPARE[2], colnames(.)))

  message("Successfully completed step 1/3...")

  ########## PART 2
  # prepare data with CNV for more convenient downstream analyses, like survival

  # First, collect events for the group 1
  # this just extracts metadata for each CNV from lesions file
  REGIONS = as.data.frame(lesions[,1:5]) %>%
    dplyr::mutate(Region = ifelse(grepl("Amplification", `Unique Name`), paste0(Descriptor, "_Amp"), paste0(Descriptor, "_Del")), .after = Descriptor)

  # Collect events for samples in group 1. Since same band can contain more than 1 peak, keep peak limits for future unique differentiation
  survival_object = 0
  # here, i is a CNV feature, anf j is a sample_id
  for(i in 1:length(REGIONS$Region)) {
    for(j in 1:1:ncol(GROUP1.CNV[,-ncol(GROUP1.CNV)])) {
      row = c(colnames(GROUP1.CNV[j]), REGIONS$Region[i], REGIONS$`Wide Peak Limits`[i], GROUP1.CNV[i,j])
      survival_object = rbind(survival_object, row)
    }
  }
  # tidy output for convenience
  survival_object = as.data.frame(survival_object)
  colnames(survival_object) = c("sample_id", "Region", "Wide Peak Limits", "Event")
  survival_object = survival_object[2:length(survival_object$sample),]
  survival_object = survival_object %>%
    dplyr::filter(!Event==0) %>%
    mutate(Event = ifelse(grepl("Amp", Region), "GAIN", "LOSS"))

  # Create final object that will be concatenated with events for second group further down
  ALL.EVENTS = survival_object

  # Collect events for samples in group 2
  survival_object = 0

  # here, i is a CNV feature, anf j is a sample_id
  for(i in 1:length(REGIONS$Region)) {
    for(j in 1:ncol(GROUP2.CNV[,-ncol(GROUP2.CNV)])) {
      row = c(colnames(GROUP2.CNV[j]), REGIONS$Region[i], REGIONS$`Wide Peak Limits`[i], GROUP2.CNV[i,j])
      survival_object = rbind(survival_object, row)
    }
  }
  # tidy output for convenience
  survival_object = as.data.frame(survival_object)
  colnames(survival_object) = c("sample_id", "Region", "Wide Peak Limits", "Event")
  survival_object = survival_object[2:length(survival_object$sample),]
  survival_object = survival_object %>% dplyr::filter(!Event==0) %>%
    mutate(Event = ifelse(grepl("Amp", Region), "GAIN", "LOSS"))

  # Now, merge with the object prepared earlier and save it as part of the output
  CNV.EVENTS = tidyr::unnest(rbind(ALL.EVENTS, survival_object))

  message("Successfully completed step 2/3...")

  ########## PART 3
  # visualize forest plot in a convenient way

  # calculate SE
  mergedPassed = GROUP1_vs_GROUP2.PASSED %>%
    dplyr::mutate(SE = (HighConfInt - LowConfInt) / 5.95)

  # order in decreasing order for better visualization
  mergedPassed = mergedPassed[order(mergedPassed$OddsRatio, decreasing = TRUE),]

  # prepare table that will be plotted
  study_table = data.frame(
    name = mergedPassed[, "Region"],
    Events_GROUP1 = paste(mergedPassed$Mutated_GROUP1, "/", mergedPassed$Mutated_GROUP1+mergedPassed$Not_Mutated_GROUP1, sep = ""),
    Events_GROUP2 = paste(mergedPassed$Mutated_GROUP2, "/", mergedPassed$Mutated_GROUP2+mergedPassed$Not_Mutated_GROUP2, sep = ""))

  # rename columns to match with comparison groups
  study_table = study_table %>%
    `colnames<-`(gsub("GROUP1", GROUPS.TO.COMPARE[1], colnames(.))) %>%
    `colnames<-`(gsub("GROUP2", GROUPS.TO.COMPARE[2], colnames(.)))

  # actual plot
  GRAPH = metaviz::viz_forest(x = mergedPassed[, c("OddsRatio", "SE")],
                                  variant = "thick",
                                  col = "Greys",
                                  xlab = "Log(OddsRatio)",
                                  annotate_CI = T,
                                  type = "study_only",
                                  study_table = study_table,
                                  text_size = text_size,
                                  table_headers = c("Region"))

  message("Successfully completed step 3/3...")

  OUTPUT = list("DISTINCT" = DISTINCT, "CNV.EVENTS" = CNV.EVENTS, "GRAPH" = GRAPH)
  return(OUTPUT)
  message("Done!")
}

#' @title Genome To Exome.
#'
#' @description Subset maf file to only features that would be available in the WEX data.
#'
#' @details To subset an incoming MAF data frame to only show features that would be available in WEX data this function was developed.
#' Pass the incoming MAF (genome) to the `maf` parameter as the only required parameter to run this function. Other parameters such as `custom_bed`,
#' `genome_build`, `padding`, and `chr_prefixed` are also available for greater control of how this function operates.
#' Refer to parameter descriptions for more information on how to use the available parameters.
#'
#' @param maf Incoming maf object. Can be maf-like data frame or maftools maf object. Required parameter. Minimum columns that should be present are Chromosome, Start_Position, and End_Position.
#' @param custom_bed Optional argument specifying a path to custom bed file for covered regions. Must be bed-like and contain chrom, start, and end position information in the first 3 columns. Other columns are disregarded if provided.
#' @param genome_build String indicating genome build of the maf file. Default is grch37, but can accept modifications of both grch37- and hg38-based builds.
#' @param padding Numeric value that will be used to pad probes in WEX data from both ends. Default is 100. After padding, overlapping features are squished together.
#' @param chr_prefixed Is the data chr-prefixed or not? Default is FALSE.
#'
#' @return A data frame of a maf-like object with the same columns as in input, but where rows are only kept for features that would be present as if the sample is WEX.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import stringr dplyr
#' @export
#'
#' @examples
#' #get all ssm in the MYC aSHM region
#' myc_ashm_maf = get_ssm_by_region(region = "8:128748352-128749427")
#' 
#' #get mutations with 100 bp padding (default)
#' maf = genome_to_exome(maf = myc_ashm_maf)
#' 
#' #get mutations covered in WEX with no padding
#' maf = genome_to_exome(maf = myc_ashm_maf, 
#'                 padding = 0) 
#' 
genome_to_exome = function(maf,
                           custom_bed,
                           genome_build = "grch37",
                           padding = 100,
                           chr_prefixed = FALSE){

  if(missing(custom_bed)){
      # first check that the genome build provided is supported
      if(! genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37", "hg38", "GRCh38", "grch38")){
        stop("The genome build specified is not currently supported. Please refer to genome build in one of the following cordinates: hg19, grch37, hs37d5, GRCh37, hg38, grch38, or GRCh38.")
      }else if(genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
        this_genome_coordinates = target_regions_grch37 # if the genome build is a flavour of hg19, get its exome space
      }else if(genome_build %in% c("hg38", "GRCh38", "grch38")){
        this_genome_coordinates = target_regions_hg38 # exome space for the variations of hg38
      }
  }else{
      this_genome_coordinates = fread(custom_bed)
      this_genome_coordinates = this_genome_coordinates[,1:3]
      colnames(this_genome_coordinates) = c("chrom", "start", "end")
  }

  # pad the ends of the probes with the padding length
  this_genome_coordinates = this_genome_coordinates %>%
    dplyr::mutate(chrom = as.character(chrom), start = start - padding, end = end + padding)

  # collapse regions if the padding results in 2 features that are overlapping
  this_genome_coordinates = this_genome_coordinates %>%
    dplyr::arrange(chrom, start) %>%
    group_by(chrom) %>%
    dplyr::mutate(indx = c(0, cumsum(as.numeric(lead(start)) > cummax(as.numeric(end)))[-n()])) %>%
    group_by(chrom, indx) %>%
    dplyr::summarise(start = first(start), end = last(end)) %>%
    dplyr::select(-indx) %>%
    ungroup %>%
    dplyr::mutate_if(is.factor, as.character)

  # handle the chr prefixes
  if(chr_prefixed & ! str_detect(this_genome_coordinates$chrom[1], "chr")){
    this_genome_coordinates = this_genome_coordinates %>%
      dplyr::mutate(chrom = paste0("chr", chrom))
  }else if(!chr_prefixed & str_detect(this_genome_coordinates$chrom[1], "chr")){
    this_genome_coordinates = this_genome_coordinates %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom, ignore.case = TRUE))
  }
  # make the coordinates as data table and set keys
  this_genome_coordinates = data.table::as.data.table(this_genome_coordinates)
  setkey(this_genome_coordinates)

  # now the maf file
  # in case the maf file is from maftools, get the ssm as data table and set keys
  if (class(maf)[1] == "MAF"){
    maf = maf@data
  }
  maf = as.data.table(maf)
  setkey(maf, Chromosome, Start_Position, End_Position)

  # subset to only features covered in exome with provided padding
  features_in_exome = foverlaps(maf, this_genome_coordinates, mult = "first", nomatch = 0L) %>% # nomatch automatically drops those that do not overlap between DTs
    as.data.frame() %>%
    dplyr::select(colnames(maf)) # make sure columns and their order is consitent with the input maf
  return(features_in_exome)
}


#' @title Tidy Lymphgen.
#'
#' @description Consolidate a column of LymphGen data in the original Subtype.Prediction output format to the GAMBLR tidy format.
#'
#' @details This function takes an incoming data frame (`df`) and consolidates a column of LymphGen data.
#' Specify the column with the lymphgen data to be processed with `lymphgen_column_in` and
#' what column to write the tidied data to with `lymphgen_column_out`.
#' In addition, the user can also run this function with `relevel = TRUE` (default is FALSE),
#' to return the output column as a factor with plot friendly levels.
#'
#' @param df Input data frame.
#' @param lymphgen_column_in The name of the column with lymphgen data to be processed.
#' @param lymphgen_column_out The name of the column to write the tidied results (optional).
#' @param relevel If TRUE, will return the output column as a factor with plot-friendly levels.
#'
#' @return A data frame with a tidied lymphGen column
#'
#' @import dplyr purrr readr stringr
#' @export
#'
#' @examples
#' metadata = get_gambl_metadata()
#' lymphgen = tidy_lymphgen(df = metadata,
#'                          lymphgen_column_in = "lymphgen_with_cnv",
#'                          lymphgen_column_out = "lymphgen_with_cnv_tidy",
#'                          relevel = TRUE)
#'
tidy_lymphgen = function(df,
                         lymphgen_column_in = "Subtype.Prediction",
                         lymphgen_column_out = "Subtype.Prediction",
                         relevel = FALSE){

  df = mutate(df, {{ lymphgen_column_out }} := case_when(
    !str_detect(.data[[lymphgen_column_in]],"/")~.data[[lymphgen_column_in]],
    str_detect(.data[[lymphgen_column_in]],"EZB")~"EZB-COMP",
    str_detect(.data[[lymphgen_column_in]],"MCD")~"MCD-COMP",
    str_detect(.data[[lymphgen_column_in]],"N1")~"N1-COMP",
    str_detect(.data[[lymphgen_column_in]],"BN2")~"BN2-COMP",
    str_detect(.data[[lymphgen_column_in]],"ST2")~"ST2-COMP"
  ))
  if(relevel){
    df <- df %>%
      mutate({
        {
          lymphgen_column_out
        }
      } := factor(
        .data[[lymphgen_column_out]],
        levels = c(
          "EZB",
          "EZB-COMP",
          "ST2",
          "ST2-COMP",
          "BN2",
          "BN2-COMP",
          "N1",
          "N1-COMP",
          "MCD",
          "MCD-COMP",
          "A53",
          "Other"
        )
      ))
  }
  return(df)
}


#' @title Consolidate Lymphgen.
#'
#' @description Replace the lymphgen column in the incoming metadata with classification for additional samples.
#'
#' @details Supplement the "lymphgen" column of the metadata with classification for additional samples.
#' Expects at least to have columns "patient_id" to bind on, and "lymphgen" to supplement the data on.
#'
#' @param sample_table Input data frame with metadata.
#' @param derived_data_path Optional argument specifying the path to a folder with files following the pattern *lymphgen.txt.
#' @param verbose Default is TRUE.
#'
#' @return A data frame with a supplemented lymphGen column.
#'
#' @import dplyr purrr readr
#' @export
#'
#' @examples
#' metadata = get_gambl_metadata()
#' consolidate_lymphgen(sample_table = metadata)
#'
consolidate_lymphgen = function(sample_table,
                                derived_data_path = "",
                                verbose = TRUE){

  if (derived_data_path == "") {
    path_to_files = check_config_value(config::get("derived_and_curated"))
    project_base = check_config_value(config::get("project_base"))
    derived_data_path = paste0(project_base, path_to_files)
    if (verbose) {
      message(
        paste0(
          "No external data path was provided, using default path ",
          derived_data_path
        )
      )
    }
  } else{
    if (verbose) {
      message(paste0("Using the specified path ", derived_data_path))
    }
  }

  lymphgen_files = dir(derived_data_path, pattern = "*lymphgen.txt")
  if (length(lymphgen_files) > 0) {
    if (verbose) {
      message(paste0(
        "Found these file(s) with lymphgen information: ",
        lymphgen_files
      ))
    }
  } else{
    if (verbose) {
      message(
        paste0(
          "No file(s) with lymphgen information were found at the path",
          lymphgen_files
        )
      )
      message(
        "If you expected the data to be present at the specified location, please ensure they follow naming convention *.lymphgen.txt"
      )
    }
  }

  for (f in lymphgen_files) {
    full = paste0(project_base, path_to_files, f)
    this_data = suppressMessages(read_tsv(full, comment = "#"))
    sample_table =
      sample_table %>% left_join(this_data, by = 'patient_id', suffix = c(".X", ".Y")) %>%
      split.default(gsub('.[XY]', '', names(.))) %>%
      map_dfc(~ if (ncol(.x) == 1)
        .x
        else
          mutate(.x, !!sym(gsub(
            '.X', '', names(.x)[1]
          )) := coalesce(!!!syms(names(
            .x
          ))))) %>%
      select(colnames(sample_table))
    sample_table = tidy_lymphgen(sample_table,
                                 lymphgen_column_in = "lymphgen",
                                 lymphgen_column_out = "lymphgen",
                                 relevel=TRUE)
  }

  return(sample_table)
}


#' @title Collate Lymphgen.
#'
#' @description Expand a sample_table (metadata) horizontally with different flavours of lymphgen data.
#'
#' @details This function takes a sample table (metadata) and adds different flavours of lymphgen data.
#' It is possible to call this function with an already subset metadata table (with sample IDs of interest) with `these_samples_metadata`.
#' If this is done, the function will join the lymphgen data with this table. Currently, the only supported `lymphgen_version` is "default".
#' For more information refer to the function examples.
#'
#' @param these_samples_metadata Optional parameter with metadata filtered for sample_ids of interest. If provided, this function will join lymphgen with this metadata, regardless of tidy TRUE/FALSE.
#' @param lymphgen_version Version of selected lymphgen, default is "default".
#' @param tidy Boolean parameter, set to TRUE for tidy format (i.e long format with no columns dropped). Default is FALSE, which returns the data in a wide format, keeping both the original Subtype. Prediction and tidied LymphGen values and puts the values from each "flavour" in its own column.
#'
#' @return A df with lymphgen information.
#'
#' @import dplyr tidyr readr stringr
#' @export
#'
#' @examples
#' library(dplyr)
#' 
#' this_meta = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "DLBCL")
#' 
#' wide_lymphgen = collate_lymphgen(these_samples_metadata = this_meta,
#'                                  lymphgen_version = "default",
#'                                  tidy = FALSE)
#'
collate_lymphgen = function(these_samples_metadata,
                            lymphgen_version = "default",
                            tidy = FALSE){

  #TODO Update the key in the config to match the version once updated, as discussed on PR.
  if(lymphgen_version == "default"){
    lymphgen_template = check_config_value(config::get("results_versioned")$lymphgen_template$default)
  }else{
    stop("Currently, only lymphgen_version = default is accepted")
  }

  #repo base
  repo_base = check_config_value(config::get("repo_base"))
  flavours = check_config_value(config::get("results_merged_wildcards")$lymphgen_template)
  flavour = str_split(flavours, pattern = ",")
  flavour = unlist(flavour)
  lymphgen_path = paste0(repo_base, lymphgen_template)

  load_lymphgen = function(flavour, lymphgen_path){
    lg_path = glue::glue(lymphgen_path)
    if(!file.exists(lg_path)){ #ignore missing flavours.
      return()
    }
    lg_df = suppressMessages(read_tsv(lg_path)) %>%
      mutate(flavour = flavour) #append the flavour in its own column called "flavour".
    return(lg_df)
  }

  lymphgen_results = lapply(flavour, load_lymphgen, lymphgen_path = lymphgen_path)
  lymphgen_results = bind_rows(lymphgen_results) #get lymphgen results tables stacked on top of each other, with the results from each flavour identified by the `flavour` column.
  lymphgen_results = tidy_lymphgen(lymphgen_results, lymphgen_column_in = "Subtype.Prediction", lymphgen_column_out = "LymphGen")
  colnames(lymphgen_results)[1] = "sample_id"

  if(!tidy){
    lymphgen_untidy = lymphgen_results %>%
      select(sample_id, Subtype.Prediction, LymphGen, flavour) %>%
      pivot_wider(names_from = flavour,
                  values_from = c(Subtype.Prediction, LymphGen),
                  names_glue = "{.value}_{flavour}")

      if(!missing(these_samples_metadata)){
        meta_data = these_samples_metadata
        lymphgen_untidy = left_join(meta_data, lymphgen_untidy)
      }

      return(lymphgen_untidy)

    }else{
    if(!missing(these_samples_metadata)){
      lymphgen_results = left_join(these_samples_metadata, lymphgen_results)
    }

    return(lymphgen_results)
  }
}


#' @title Collate Quality Control Results.
#'
#' @description Expand a metadata table horizontally with quality control metrics.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table df with sample ids in the first column.
#' @param seq_type_filter default is genome, capture is also available for unix_group icgc_dart.
#'
#' @return The sample table with additional columns.
#'
#' @import dplyr readr
#' 
#' @noRd
#'
#' @examples
#' qc_metrics = collate_qc_results(sample_table = sample_table)
#'
collate_qc_results = function(sample_table,
                              seq_type_filter = "genome"){

  if(! seq_type_filter %in% c("genome", "capture")){
    stop("Please provide a valid seq_type (\"genome\" or \"capture\").")
  }

  #get paths
  base = check_config_value(config::get("project_base"))
  qc_template = check_config_value(config::get("qc_met"))

  #icgc_dart
  unix_group = "icgc_dart"
  icgc_qc_path = glue::glue(qc_template)
  icgc_qc_path_full = paste0(base, icgc_qc_path)

  #gambl
  unix_group = "gambl"
  gambl_qc_path = glue::glue(qc_template)
  gambl_qc_path_full = paste0(base, gambl_qc_path)

  #read in icgc qc data, rename sample id column and filter on samples in sample id in sample_table
  icgc_qc = suppressMessages(read_tsv(icgc_qc_path_full)) %>%
      dplyr::rename(sample_id = UID) %>%
      dplyr::select(-SeqType)

  #read in gambl qc data (if seq_type_filter set to "genome"), rename sample id column and filter on samples in sample id in sample_table
  if(seq_type_filter == "genome"){
    gambl_qc = suppressMessages(read_tsv(gambl_qc_path_full)) %>%
      dplyr::rename(sample_id = UID) %>%
      dplyr::select(-SeqType)

    #join gambl and icgc QC data
    full_qc = rbind(gambl_qc, icgc_qc)
    sample_table = left_join(sample_table, full_qc)

    #print n samples with QC metrics
    qc_samples = length(unique(full_qc$sample_id))
    message(paste("QC metrics for", qc_samples, "samples retrieved."))

  }else{
    message("Currently, seq_type_filter = \"capture\" is only available for unix_group \"icgc_dart\". Only QC metrics for icgc_dart will be returned.")
    #TO DO: Remove this once capture metrics have been generated for gambl samples.
    sample_table = left_join(sample_table, icgc_qc)

    #print n samples with QC metrics
    qc_samples = length(unique(icgc_qc$sample_id))
    message(paste("QC metrics for", qc_samples, "samples retrieved."))
  }
  return(sample_table)
}


#' @title Standardize Chromosome Prefix.
#'
#' @description Standardize the chr prefix in a vector of chromosome names based on projection.
#'
#' @details INTERNAL FUNCTION, not meant for out-of-package use.
#'
#' @param incoming_vector Input vector of any length with chromosome names.
#' @param projection Projection to which chr prefix should be standardized.
#'
#' @return A vector of chromosome names with prefix standardized to projection
#' 
#' @noRd
#'
#' @examples
#' these_chrs = c(8, "13", "chr4", "chrY")
#' 
#' standardize_chr_prefix(incoming_vector = these_chrs,
#'                        projection = "hg38")
#'
standardize_chr_prefix = function(incoming_vector,
                                  projection){

  if (projection %in% c("grch37", "grch38")) {
    output_vector = gsub("chr", "", incoming_vector)
  } else {
    output_vector = gsub("chr", "", incoming_vector) # if there is a mix of prefixed and non-prefixed options
    output_vector = paste0("chr", output_vector)
  }
  return(output_vector)
}


#' @title Calculate proportion of genome altered by CNV.
#'
#' @description `calculate_pga` returns a data.frame with estimated proportion of genome altered for each sample.
#'
#' @details This function calculates the percent of genome altered (PGA) by CNV. It takes into account the total length of
#' sample's CNV and relates it to the total genome length to return the proportion affected by CNV. The input is expected to be a seg file.
#' The path to a local SEG file can be provided instead. If The custom seg file is provided, the minimum required columns are
#' sample, chrom, start, end, and log.ratio. The function can work with either individual or multi-sample seg files. The telomeres are always
#' excluded from calculation, and centromeres/sex chromosomes can be optionally included or excluded.
#'
#' @param this_seg Input data frame of seg file.
#' @param seg_path Optionally, specify the path to a local seg file.
#' @param projection Argument specifying the projection of seg file, which will determine chr prefix, chromosome coordinates, and genome size. Default is grch37, but hg38 is also accepted.
#' @param cutoff The minimum log.ratio for the segment to be considered as CNV. Default is 0.56, which is 1 copy. This value is expected to be a positive float of log.ratio for both deletions and amplifications.
#' @param exclude_sex Boolean argument specifying whether to exclude sex chromosomes from calculation. Default is TRUE.
#' @param exclude_centromeres Boolean argument specifying whether to exclude centromeres from calculation. Default is TRUE.
#'
#' @return data frame
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr readr tidyr
#' @export
#'
#' @examples
#' library(dplyr)
#' 
#' sample_seg = get_sample_cn_segments(this_sample_id = "14-36022T") %>%
#'  rename("sample"="ID")
#' 
#' calculate_pga(this_seg = sample_seg)
#' 
#' calculate_pga(this_seg = sample_seg,
#'               exclude_sex = FALSE)
#'
#' multi_sample_seg = rbind(get_sample_cn_segments(this_sample_id = "14-36022T"),
#'                          get_sample_cn_segments(this_sample_id = "BLGSP-71-21-00243-01A-11E")) %>%
#'                          rename("sample"="ID")
#' 
#' calculate_pga(this_seg = multi_sample_seg)
#'
calculate_pga = function(this_seg,
                         seg_path,
                         projection = "grch37",
                         cutoff = 0.56,
                         exclude_sex = TRUE,
                         exclude_centromeres = TRUE) {
  # check for required argument
  if (missing(this_seg) & missing (seg_path)) {
    stop("Please provide the data frame of seg file or path to the local seg.")
  }

  # ensure the specified projection is correct and define chromosome coordinates
  if (projection == "grch37") {
    chr_coordinates = chromosome_arms_grch37
  } else if (projection == "hg38") {
    chr_coordinates = chromosome_arms_hg38
  } else {
    stop(
      "You specified projection that is currently not supported. Please provide seg files in either hg38 or grch37."
    )
  }

  # exclude sex chromosomes
  if (exclude_sex) {
    chr_coordinates = chr_coordinates %>%
      dplyr::filter(!grepl("X|Y", chromosome))
  }

  # does the user's seg file contain centromeres?
  if (exclude_centromeres) {
    chr_coordinates = chr_coordinates %>%
      group_by(chromosome) %>%
      mutate(start = min(start),
             end = max(end)) %>%
      ungroup
  }

  # total size of genome in this projection
  genome_size = chr_coordinates %>%
    mutate(size = end - start) %>%
    summarise(genome_size = sum(size)) %>%
    pull(genome_size)

  # prepare for the overlaps
  chr_coordinates = as.data.table(chr_coordinates)  %>%
    rename("arm_start" = "start",
           "arm_end" = "end")
  setkey(chr_coordinates, chromosome, arm_start, arm_end)

  # work out the seg file
  if (!missing(seg_path)) {
    message(paste0("Reading thhe seg file from ", seg_path))
    this_seg = suppressMessages(read_tsv(seg_path))
  }

  # preserve the sample ids to account later for those with 0 PGA
  sample_set = this_seg %>% distinct(sample)

  this_seg = this_seg %>%
    dplyr::filter(abs(log.ratio) >= cutoff) %>%
    dplyr::relocate(sample, .after = last_col())

  # ensure consistent chromosome prefixing
  if (projection == "grch37") {
    this_seg$chrom = gsub("chr", "", this_seg$chrom)
  } else {
    this_seg$chrom = gsub("chr", "", this_seg$chrom) # if there is a mish-mash of prefixes, strip them all
    this_seg$chrom = paste0("chr", this_seg$chrom)
  }

  # exclude sex chromosomes
  if (exclude_sex) {
    this_seg = this_seg %>%
      dplyr::filter(!grepl("X|Y", chrom))
  }

  # prepare for the overlaps
  this_seg = as.data.table(this_seg)
  setkey(this_seg, chrom, start, end)

  # what are the segments that overlap good regions in chromosome coordinates?
  this_seg = foverlaps(
    this_seg,
    chr_coordinates,
    by.x = c("chrom", "start", "end"),
    by.y = c('chromosome', 'arm_start', 'arm_end'),
    nomatch = 0L
  ) %>%
    as.data.frame %>%
    arrange(sample, chrom, start)

  # calculate total length of CNV
  affected_regions = this_seg %>%
    dplyr::mutate(size = end - start) %>%
    group_by(sample) %>%
    summarise(total = sum(size))

  affected_regions$PGA = affected_regions$total / genome_size

  # now add any samples that can have 0 PGA
  affected_regions = base::merge(sample_set,
                                 affected_regions,
                                 all.x = TRUE) %>%
    replace_na(list(PGA = 0))

  affected_regions = affected_regions %>%
    select(-total) %>%
    `names<-`(c("sample_id", "PGA")) %>%
    as.data.frame()

  return(affected_regions)

}


#' @title Adjust ploidy for samples with CNV data.
#'
#' @description `adjust_ploidy` returns a seg file with log.ratios adjusted to the overall sample ploidy.
#'
#' @details This function adjusts the ploidy of the sample using the percent of genome altered (PGA). The PGA is calculated internally, but can also be optionally provided as data frame
#' if calculated from other sources. Only the samples above the threshold-provided PGA will have ploidy adjusted. The function can work with either individual or
#' multi-sample seg file. The telomeres are always excluded from calculation, and sex chromosomes can be optionally included or excluded. The supported projections are grch37 and hg38.
#' The chromosome prefix is handled internally per projection and does not need to be consistent.
#'
#' @param this_seg Input data frame of seg file.
#' @param seg_path Optionally, specify the path to a local seg file.
#' @param projection Argument specifying the projection of seg file, which will determine chr prefix and genome size. Default is grch37, but hg38 is also accepted.
#' @param pga If PGA is calculated through other sources, the data frame with columns sample_id and PGA can be provided in this argument.
#' @param pga_cutoff Minimum PGA for the sample to adjust ploidy. Default is 0.05 (5%).
#' @param exclude_sex Boolean argument specifying whether to exclude sex chromosomes from calculation. Default is TRUE.
#' @param return_seg Boolean argument specifying whether to return a data frame in seg-consistent format, or a raw data frame with all step-by-step transformations. Default is TRUE.
#'
#' @return A data frame in seg-consistent format with ploidy-adjusted log ratios.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr tidyr readr
#' @export
#'
#' @examples
#' library(dplyr)
#' 
#' sample_seg = get_sample_cn_segments(this_sample_id = "14-36022T") %>%
#'  rename("sample"="ID")
#' 
#' adjust_ploidy(this_seg = sample_seg)
#'
#' multi_sample_seg = rbind(get_sample_cn_segments(this_sample_id = "14-36022T"),
#'                          get_sample_cn_segments(this_sample_id = "BLGSP-71-21-00243-01A-11E")) %>%
#'                          rename("sample" = "ID")
#' 
#' adjust_ploidy(this_seg = multi_sample_seg)
#'
adjust_ploidy = function(this_seg,
                         seg_path,
                         projection = "grch37",
                         pga,
                         pga_cutoff = 0.05,
                         exclude_sex = TRUE,
                         return_seg = TRUE) {

  # ensure the specified projection is correct
  # this is only needed to adjust ploidy if the pre-adjusted ploidy is not provided
  if (!projection %in% c("grch37", "hg38") & missing(pga)) {
    stop(
      "You specified projection that is currently not supported. Please provide seg files in either hg38 or grch37."
    )
  }

  # if the seg is a local file, read it in
  if (!missing(seg_path)) {
    message(paste0("Reading thhe seg file from ", seg_path))
    this_seg = suppressMessages(read_tsv(seg_path))
  }

  # ensure consistent chromosome prefixing
  if (projection == "grch37") {
    this_seg$chrom = gsub("chr", "", this_seg$chrom)
  } else {
    this_seg$chrom = gsub("chr", "", this_seg$chrom) # if there is a mish-mash of prefixes, strip them all
    this_seg$chrom = paste0("chr", this_seg$chrom)
  }

  # exclude sex chromosomes
  if (exclude_sex) {
    this_seg = this_seg %>%
      dplyr::filter(!grepl("X|Y", chrom))
  }

  # if PGA is called with custom parameters or obtained elsewhere, use it - but if not, calculate on-the-fly
  if (missing(pga)) {
    message("Calculating PGA ...")
    pga = calculate_pga(this_seg = this_seg,
                        projection = projection,
                        exclude_sex = exclude_sex)
  } else {
    # ensure the column is named "sample_id" if it came outside GAMBLR
    if (!"sample_id" %in% colnames(pga)) {
      stop("Please ensure the column with sample ids in your PGA data frame is named sample_id.")
    }
  }

  # add the PGA information to the seg file
  this_seg = left_join(this_seg,
                       pga,
                       by = c("sample" = "sample_id"))

  # By how much we should adjust the ploidy?
  this_seg = this_seg %>%
    group_by(sample) %>% # account for multi-sample seg files
    dplyr::mutate(
      adjust = mean((2 * 2 ^ log.ratio)),
      # convert log.ratio to absolute CN and find average
      adjust = PGA * adjust,
      # what is the ploidy of genome affected by CNV?
      neutral = (1 - PGA) * 2,
      # how much of the genome is not affected by CNV? Assume it's in diploid state
      adjust = adjust + neutral,
      # overall genome's ploidy status
      adjust = abs(adjust - 2)
    ) %>% # by how much the ploidy should be adjusted? From average sample ploidy take out the diploid state
    dplyr::mutate(need_to_adjust = ifelse((PGA > pga_cutoff &
                                             !log.ratio == 0), "TRUE", "FALSE")) # We only need to adjust ploidy if the PGA is above cut-off and segment is not diploid

  # Adjust ploidy
  this_seg = this_seg %>%
    ungroup() %>%
    mutate(cn = (2*2^log.ratio), .before = log.ratio) %>% # convert log.ratios to absolute CN
    mutate(new_cn = ifelse(need_to_adjust == "TRUE", # always round to the integer, but adjust only if needed
                           round(cn-adjust),
                           round(cn))) %>%
    mutate(new_log.ratio = ifelse(need_to_adjust == "TRUE", # log transform the new CN states
                                  (log(abs(new_cn),2)-1),
                                  log.ratio)) %>%
    mutate(new_log.ratio = ifelse(new_log.ratio == -Inf,-10, new_log.ratio)) # deal with the -Inf for low negative numbers

  # this allows user to see all the transformations if they wish or return the standard seg file
  if (return_seg) {
    message("Returning the seg file with ploidy-adjusted CN ...")
    this_seg = this_seg %>%
      mutate(log.ratio = new_log.ratio) %>% # assign the new log.ratio
      select(sample, chrom, start, end, LOH_flag, log.ratio) # select only columns of the seg file
  }


  return(this_seg)
}


#' @title Subset CN States.
#'
#' @description Get the available CN states in the incoming data frame.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::fancy_multisample_ideo], for sub-setting copy number information based on segments available in cn data
#'
#' @param cn_segments DF with copy number segments, usually retrieved from get_sample_cn_segments.
#' @param include_2 Optional parameter for including or omit CN state == 2. Default is FALSE.
#' @param samplen Numeric value that annotates the sample order.
#'
#' @return Nothing.
#' 
#' @noRd
#'
#' @examples
#' cn_states = get_sample_cn_segments(multiple_samples = TRUE,
#'                                    sample_list = c("00-15201_tumorA",
#'                                                    "HTMCP-01-06-00422-01A-01D"),
#'                                    streamlined = FALSE)
#' 
#' subset_cnstates(cn_segments = cn_states,
#'                 samplen = 1)
#'
subset_cnstates = function(cn_segments,
                           include_2 = FALSE,
                           samplen){

  #transform CN states > 6 = 6 (to reflect the current copy number palette in gamblr)
  cn_segments$CN[cn_segments$CN > 6] = 6

  #filter out CN == 2
  if(!include_2){
    cn_segments = subset(cn_segments, CN != 2)
  }

  #update CN annotations (if present in cn_segment data).
  cn_segments$CN = paste0("cn_", cn_segments$CN , "_sample", samplen)

  #convert to factor.
  cn_segments$CN = as.factor(cn_segments$CN)

  #split cn_segments on available factors and lists into the global environment.
  l = split(cn_segments, cn_segments$CN)
  list2env(l, envir = .GlobalEnv)
}


#' @title Compare segmented data for multiple samples.
#'
#' @description `cnvKompare` returns a list in variable data formats allowing to evaluate concordance of CNV data between multiple samples.
#'
#' @details This function will compare CNV data between samples with multiple time points. It can also handle same-sample comparison
#' between different CNV callers if the sample ID is specified in unique fashion. For groups with more than 2 samples,
#' optionally pairwise comparisons can be performed. The comparison is made based on the internally calculated score,
#' which reflects the percentage of each cytoband covered by CNV (rounded to the nearest 5%) and its absolute CN. Optionally,
#' the heatmap of cnvKompare scores can be returned. In addition, the function will return all concordant and discordant cytobands.
#' Finally, the time series plot of CNV log ratios will be returned for all lymphoma genes, with further functionality to subset
#' it to a panel of genes of interest.
#'
#' @param patient_id Specify patient_id to retrieve sample ids from GAMBL metadata.
#' @param these_sample_ids Optionally, specify sample ids for comparison.
#' @param this_seg Optional input data frame of seg file. Must adhere to seg format.
#' @param seg_path Optionally, specify the path to a local seg file. Must adhere to seg format.
#' @param genes_of_interest Provide specific genes to be displayed on the time-series plot.
#' @param projection Argument specifying the projection of seg file, which will determine coordinates of the cytobands. Default is grch37, but hg38 is also accepted.
#' @param ignore_cytoband_labels Cytobands to be ignored. By default, "acen", "gvar", "stalk" are excluded.
#' @param max_overlap For a time-series plot, how many maximum overlapping points are allowed?
#' @param min_concordance Integer value from 0 to 100 to indicate the minimum required similarity between cytobands to be considered concordant. The default is 90 (90%).
#' @param exclude_sex Boolean argument specifying whether to exclude sex chromosomes from calculation. Default is FALSE.
#' @param return_heatmap Boolean argument specifying whether to return a heatmap of cnvKompare scores. Default is TRUE.
#' @param compare_pairwise Boolean argument specifying whether to perform pairwise comparisons if there are more than 2 time points in the group. Default is TRUE.
#'
#' @return A list of overall and pairwise percent concordance, concordant and discordant cytobands, comparison heatmap of cnvKompare scores, and time series ggplot object.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr tidyr circlize ComplexHeatmap ggplot2 ggrepel readr tibble
#' @importFrom plyr round_any
#' @export
#'
#' @examples
#' cnvKompare(patient_id = "13-26835",
#'            genes_of_interest = c("EZH2",
#'                                  "TP53",
#'                                  "MYC",
#'                                  "CREBBP",
#'                                  "GNA13"),
#'            projection = "hg38")
#' 
cnvKompare = function(patient_id,
                      these_sample_ids,
                      this_seg,
                      seg_path,
                      genes_of_interest,
                      projection = "grch37",
                      ignore_cytoband_labels = c("acen", "gvar", "stalk"),
                      max_overlap = 20,
                      min_concordance = 90,
                      exclude_sex = FALSE,
                      return_heatmap = TRUE,
                      compare_pairwise = TRUE) {
  # initialize output list
  output = list()

  # convert min concordance pct to fraction
  min_concordance = min_concordance/100

  # check that sample identifiers are provided
  if (missing(patient_id) & missing(these_sample_ids)) {
    stop("Please provide patient id or sample ids for comparison.")
  }

  # retrieve sample ids if only patient id is specified
  if (missing(these_sample_ids)) {
    these_sample_ids = get_gambl_metadata()
    these_sample_ids = dplyr::filter(these_sample_ids, patient_id == {{ patient_id }})
    these_sample_ids = pull(these_sample_ids, sample_id)
    message(paste0(
      "Found ",
      length(these_sample_ids),
      " samples for patient ",
      patient_id,
      " ..."
    ))
  }

  # get cytobands
  if (projection %in% c("hg19", "grch37")) {
    cytobands = circlize::read.cytoband(species = "hg19")$df %>%
      mutate(V1 = gsub("chr", "", V1))
  } else if (projection %in% c("hg38", "grch38")) {
    cytobands = circlize::read.cytoband(species = "hg38")$df
  } else {
    stop("Please specify one of hg19, grch37, hg38, or grch38 projections.")
  }
  cytobands = cytobands %>%
    `names<-`(c("cb.chromosome", "cb.start", "cb.end", "cb.name", "label")) %>%
    dplyr::filter(!label %in% ignore_cytoband_labels)
  if (exclude_sex) {
    cytobands = dplyr::filter(cytobands,!grepl("X|Y", cb.chromosome))
  }
  cytobands = as.data.table(cytobands)
  setkey(cytobands, cb.chromosome, cb.start, cb.end)

  # get the multi-sample seg file
  if (!missing(seg_path)) {
    these_samples_seg = suppressMessages(read_tsv(seg_path)) %>%
      `names<-`(c(ID, chrom, start, end, LOH_flag, log.ratio)) %>%
      dplyr::mutate(CN = (2 * 2 ^ log.ratio))
  } else if (!missing(this_seg)) {
    these_samples_seg = this_seg %>%
      `names<-`(c(ID, chrom, start, end, LOH_flag, log.ratio)) %>%
      dplyr::mutate(CN = (2 * 2 ^ log.ratio))
  } else {
    message("Retreiving the CNV data using GAMBLR ...")
    these_samples_seg = get_sample_cn_segments(multiple_samples = TRUE,
                                               sample_list = these_sample_ids,
                                               from_flatfile = TRUE,
                                               projection = projection,
                                               with_chr_prefix = TRUE)
  }

  these_samples_seg = these_samples_seg  %>%
    dplyr::filter(ID %in% these_sample_ids) %>% # if user-provided seg, ensure only samples of comparison are present
    relocate(ID, .after = last_col())
  if (exclude_sex) {
    these_samples_seg = dplyr::filter(these_samples_seg,!grepl("X|Y", chrom))
  }

  these_samples_seg = as.data.table(these_samples_seg)
  setkey(these_samples_seg, chrom, start, end)

  # overlap seg with cytoband regions
  # if segment extends beyond cytoband, cut it at the cytoband coordinates
  cytoband_overlap =
    foverlaps(cytobands,
              these_samples_seg,
              nomatch = 0) %>%
    as.data.frame %>%
    group_by(cb.name) %>%
    mutate(
      start = ifelse(start > cb.start, start, cb.start),
      end = ifelse(end < cb.end, end, cb.end)
    ) %>%
    ungroup

  # calculate % of cytoband covered by CNV and concordance score
  message("Calculating CNV concordance ...")
  for_output =
    cytoband_overlap %>%
    mutate(
      band_length = cb.end - cb.start,
      cnv_length = end - start,
      # round the % of cytoband covered to the nearest 5%
      pct_covered = plyr::round_any((cnv_length / band_length * 100), 5, f = round)
    ) %>%
    # name each cytoband as chr_cytoband
    mutate(name = paste0(cb.chromosome, "_", cb.name)) %>%
    # calculate % covered by cytoband and it's CN
    group_by(ID, name) %>%
    mutate(pct_covered = sum(pct_covered),
           CN = mean(CN)) %>%
    ungroup() %>%
    arrange(name, ID) %>%
    distinct(ID, CN, pct_covered, name, .keep_all = TRUE) %>%
    # designate score of %covered*CN to infer intra-sample concordance
    mutate(score = CN * pct_covered)

  # overall concordance
  concordance =
    for_output %>%
    select(ID, name, score) %>%
    spread(., ID, score) %>%
    column_to_rownames("name")

  overall_concordance = ifelse((rowSums(concordance) >= ((concordance[, 1] * ncol(concordance))*min_concordance) &
                                  rowSums(concordance) <= ((concordance[, 1] * ncol(concordance))*(1+(1-min_concordance)))),
                               "YES",
                               "NO")
  overall_concordance = overall_concordance[!is.na(overall_concordance)]
  overall_concordance_pct = round(((
    sum(overall_concordance == "YES") / length(overall_concordance)
  ) * 100), 2)
  output$overall_concordance_pct = overall_concordance_pct

  # return cytobands consistent across samples
  concordant_cytobands =
    for_output %>%
    # output-specific
    select(ID, cb.chromosome, cb.start, cb.end, name, score, log.ratio) %>%
    dplyr::filter(name %in% names(overall_concordance[overall_concordance == "YES"]))

  output$concordant_cytobands = concordant_cytobands

  # return cytobands discordant across samples
  discordant_cytobands =
    for_output %>%
    # output-specific
    select(ID,
           cb.chromosome,
           cb.start,
           cb.end,
           name,
           pct_covered,
           log.ratio,
           score) %>%
    dplyr::filter(name %in% names(overall_concordance[overall_concordance == "NO"]))

  output$discordant_cytobands = discordant_cytobands

  # heatmap of cnvKompare scores
  if (return_heatmap) {
    message("Building heatmap ...")
    hmap_legend_param = list(title = "cnvKompare score")
    hMap = concordance %>%
      as.matrix() %>%
      t %>%
      ComplexHeatmap::Heatmap(
        .,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        heatmap_legend_param = hmap_legend_param
      )
    output$Heatmap = hMap
  }

  # VAF-like plot
  # genes
  if (projection %in% c("hg19", "grch37")) {
    for_plot_lg = grch37_lymphoma_genes_bed %>%
      as.data.table()
  } else if (projection %in% c("hg38", "grch38")) {
    for_plot_lg = hg38_lymphoma_genes_bed %>%
      as.data.table()
  }
  # did user specify particular genes of interest to display on the plot?
  if (!missing(genes_of_interest)) {
    message("Subsetting lymphoma genes to specified genes of interest ...")
    for_plot_lg = dplyr::filter(for_plot_lg, hgnc_symbol %in% genes_of_interest)
  }
  setkey(for_plot_lg, chromosome_name, start_position, end_position)

  # CNV data
  for_plot = as.data.table(for_output)
  setkey(for_plot, cb.chromosome, start, end)

  # generate plot
  time_plot =
    foverlaps(for_plot_lg,
              for_plot,
              nomatch = NULL) %>%
    as.data.frame %>%
    select(ID, hgnc_symbol, log.ratio) %>%
    ggplot(., aes(x = ID, y = log.ratio, group = hgnc_symbol)) +
    geom_line(aes(color = hgnc_symbol)) +
    ggrepel::geom_label_repel(
      aes(label = hgnc_symbol, color = hgnc_symbol),
      nudge_x = 0.1,
      na.rm = TRUE,
      label.size = NA,
      fill = NA,
      segment.color = "transparent",
      max.overlaps = max_overlap
    ) +
    theme_Morons() +
    theme(legend.position = "none",
          axis.title.x = element_blank())

  output$time_plot = time_plot

  # for groups with >2 samples, make pairwise comparisons
  if (compare_pairwise & length(these_sample_ids) > 2) {
    message("Performing pairwise comparisons ...")

    # generate all possible combinations
    possible_combinations = apply(combn(these_sample_ids, 2), 2, paste, collapse =
                                    '--')

    for (combination in possible_combinations) {
      # samples in this pairwise comparison
      these_samples = unlist(strsplit(combination, split = "--"))

      # pct concordance in this pair
      this_concordance = concordance[, these_samples]
      this_concordance = ifelse((this_concordance[, 1] >= (this_concordance[, 2])*min_concordance) &
                                  this_concordance[, 1] <= (this_concordance[, 2])*(1+(1-min_concordance)),
                                "YES",
                                "NO")
      names(this_concordance) = rownames(concordance)
      this_concordance = this_concordance[!is.na(this_concordance)]
      this_concordance_pct = round(((
        sum(this_concordance == "YES") / length(this_concordance)
      ) * 100), 2)
      output$pairwise_comparisons[[combination]]$pairwise_concordance_pct = this_concordance_pct

      # return cytobands consistent in this pair
      concordant_cytobands =
        for_output %>%
        dplyr::filter(ID %in% these_samples) %>%
        # output-specific
        select(ID, cb.chromosome, cb.start, cb.end, name, score) %>%
        dplyr::filter(name %in% names(this_concordance[this_concordance == "YES"]))

      output$pairwise_comparisons[[combination]]$concordant_cytobands = concordant_cytobands

      # return cytobands discordant in this pair
      discordant_cytobands =
        for_output %>%
        dplyr::filter(ID %in% these_samples) %>%
        # output-specific
        select(ID,
               cb.chromosome,
               cb.start,
               cb.end,
               name,
               pct_covered,
               log.ratio,
               score) %>%
        dplyr::filter(name %in% names(this_concordance[this_concordance == "NO"]))

      output$pairwise_comparisons[[combination]]$discordant_cytobands = discordant_cytobands

    }
  }
  message("DONE!")
  return(output)
}


#' @title Cleanup MAF.
#'
#' @description Transform input maf columns to allow for usage of dplyr verbs.
#'
#' @details Transform input maf columns to allow for usage of dplyr verbs.
#' Allowing for a stright-forward plotting workflow as well as downstream data aggregation and manipulation.
#' This function expects a set number of columns to exist in the incoming maf in order for this to work.
#' To view the columns, see bundled data; [GAMBLR::grande_maf].
#'
#' @param maf_df input MAF data frame.
#'
#' @return maf_df with transformed columns
#'
#' @import dplyr
#' @export
#'
#' @examples
#' 
#' clean_maf = cleanup_maf(maf_df = grande_maf)
#'
cleanup_maf = function(maf_df){

  #cleanup various columns that store text to make them useful (make numeric, drop denominator etc)
  maf_df = mutate(maf_df,EXON = gsub("/.+", "", EXON)) %>%
    mutate(EXON = as.numeric(EXON)) %>%
    mutate(INTRON = gsub("/.+", "", INTRON)) %>%
    mutate(INTRON = as.numeric(INTRON)) %>%
    mutate(CDS_position = gsub("/.+", "", CDS_position)) %>%
    mutate(CDS_position = as.numeric(as.character(CDS_position))) %>%
    mutate(cDNA_position = gsub("/.+", "", cDNA_position)) %>%
    mutate(cDNA_position = as.numeric(as.character(cDNA_position))) %>%
    mutate(Protein_position = gsub("/.+", "", Protein_position)) %>%
    mutate(Protein_position = as.numeric(as.character(Protein_position)))

  return(maf_df)
}


#' @title Supplement MAF.
#'
#' @description Complement maf with missing samples.
#'
#' @details Specify the initial MAF with `incoming_maf` (to be supplemented with missing samples) and
#' give the function a filtered metadata table (with the sample IDs of interest) to the `these_samples_metadata`.
#'
#' @param incoming_maf The initial MAF data frame to be supplemented with missing samples.
#' @param these_samples_metadata The metadata data frame that contains Tumor_Sample_Barcode column with ids to be present in the complemented MAF.
#'
#' @return maf_df with complemented Tumor_Sample_Barcode and other columns ready to be used downstream.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' 
#' small_maf = get_coding_ssm(limit_cohort = "dlbcl_reddy",
#'                            seq_type = "capture") %>% 
#'  dplyr::filter(Hugo_Symbol=="MYC")
#' 
#' reddy_meta = get_gambl_metadata(seq_type_filter = "capture") %>% 
#'  dplyr::filter(cohort=="dlbcl_reddy")
#' 
#' complete_maf = supplement_maf(incoming_maf = small_maf,
#'                               these_samples_metadata = reddy_meta)
#'
supplement_maf <- function(incoming_maf,
                           these_samples_metadata){

  missing_sample_ids = setdiff(these_samples_metadata$Tumor_Sample_Barcode,
                               incoming_maf$Tumor_Sample_Barcode)

  missing_sample_maf = incoming_maf %>%
    dplyr::filter(Tumor_Sample_Barcode == "Imaginary Sample ID") %>%
    add_row(Tumor_Sample_Barcode = missing_sample_ids,
           Hugo_Symbol = "GARBAGE",
           Chromosome = ifelse(stringr::str_detect(incoming_maf$Chromosome[1], "chr"), "chr1", "1"),
           Start_Position = 1,
           End_Position = 1,
           Variant_Classification = "Missense_Mutation")

  full_maf = rbind(incoming_maf, missing_sample_maf)
  return(full_maf)
}
