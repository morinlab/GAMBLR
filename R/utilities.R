#global environment
coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")
rainfall_conv = c("T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G", "InDel")
names(rainfall_conv) = c('A>G', 'T>C', 'C>T', 'G>A', 'A>T', 'T>A', 'A>C', 'T>G', 'C>A', 'G>T', 'C>G', 'G>C', 'InDel')
ssh_session <<- NULL


#' Check if code is running remotely and (optionally) attempt a connection and set global ssh_session variable
#'
#' @param auto_connect Set to TRUE if you want the function to create an ssh session (if necessary)
#'
#' @return NULL. This function makes a new ssh session, if necessary, and stores it in a global variable (ssh_session) 
#' @export
#'
#' @examples ssh_session = get_ssh_session()
check_host = function(auto_connect=FALSE){
  hostname = Sys.info()["nodename"]
  if(grepl("bcgsc.ca",hostname)){
    #we are on the GSC network
  }else{
    # we are on some computer not on the GSC network (needs ssh_session)
    if(class(ssh_session)=="ssh_session"){
      message("active ssh session detected")
      return()
    }else{
      if(auto_connect){
        session = get_ssh_session()
        assign("ssh_session", session, envir = .GlobalEnv)
      }else{
        message("You appear to be using GAMBLR on your local computer. Be sure to set up an ssh session!")
        message("?GAMBLR::get_ssh_session for more info")
        
      }
    }
  }
}

#' Count the variants in a region with a variety of filtering options
#'
#' @param region 
#' @param chromosome 
#' @param start 
#' @param end 
#' @param these_samples_metadata 
#' @param count_by Defaults to counting all variants. Specify 'sample_id' if you want to collapse and count only one per sample
#' @param seq_type 
#' @param ssh_session 
#'
#' @return
#' @export
#'
#' @examples
count_ssm_by_region = function(region,chromosome,start,end,these_samples_metadata,count_by,seq_type="genome",ssh_session){
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata(seq_type_filter=seq_type)
  }
  if(missing(region)){
    region_muts = get_ssm_by_region(chromosome=chromosome,qstart=start,qend=end,streamlined = TRUE,ssh_session = ssh_session)
  }else{
    region_muts = get_ssm_by_region(region=region,streamlined = TRUE,ssh_session = ssh_session)
  }
  keep_muts = left_join(region_muts,these_samples_metadata) %>% dplyr::filter(!is.na(sample_id))
  if(missing(count_by)){
    #count everything even if some mutations are from the same patient
    return(nrow(keep_muts))
  }else if(count_by == "sample_id"){
    return(length(unique(keep_muts$sample_id)))
  }else{
    print("Not sure what to count")
  }
}

#' Split a contiguous genomic region on a chromosome into non-overlapping bins
#'
#' @param chromosome 
#' @param start 
#' @param end 
#' @param bin_size 
#'
#' @return Data frame describing the bins various ways
#' @export
#'
#' @examples 
#' chr8q_bins = region_to_bins(chromosome="8",start=48100000,end=146364022,bin_size = 20000)
region_to_bins = function(chromosome="chr1",start=10000,end=121500000,bin_size=2000){
  bin_df = data.frame(bin_chr = chromosome,bin_start=seq(start,end,bin_size))
  bin_df = mutate(bin_df,bin_end = bin_start+ bin_size) %>%
    dplyr::filter(bin_end<=end) %>% 
    mutate(region = paste0(bin_chr,":",bin_start,"-",bin_end))
  
  return(bin_df)
  
}

#' Create an ssh session to the GSC (requires active VPN connection)
#'
#' @param host (default is gphost01.bcgsc.ca)
#'
#' @return
#' @export
#'
#' @examples
get_ssh_session = function(host="gphost01.bcgsc.ca"){
  if (!requireNamespace("ssh", quietly = TRUE)) {
    warning("The ssh package must be installed to use this functionality")
    #Either exit or do something that does not require ssh
    return(NULL)
  }
  session = ssh::ssh_connect(host=host)
  return(session)
}


#' Get regions from genes.
#'
#' @param gene_symbol Gene symbol (e.g BCL2).
#' @param ensembl_id Esemble ID (e.g ENSG00000171791).
#' @param genome_build Reference genome build.
#' @param return_as Specify the type of return. Default is bed (chr:start-end), other acceptable inputs are "df".
#'
#' @return
#' @export
#'
#' @examples
#' bcl2_region = gene_to_region(gene_symbol = "BCL2", genome_build = "grch37")
#' bcl2_region = gene_to_region(ensembl_id = "ENSG00000171791", genome_build = "grch37")
#'
gene_to_region = function(gene_symbol,
                          ensembl_id,
                          genome_build = "grch37",
                          return_as = "bed"){

  if(genome_build == "grch37"){
    if(!missing(gene_symbol) && missing(ensembl_id)){
      gene_coordinates = dplyr::filter(grch37_gene_coordinates, hugo_symbol %in% gene_symbol)
    }

    if(missing(gene_symbol) && !missing(ensembl_id)){
      gene_coordinates = dplyr::filter(grch37_gene_coordinates, ensembl_gene_id %in% gene_symbol)
    }
  }

  if(genome_build == "hg38"){
    if(!missing(gene_symbol) && missing(ensembl_id)){
      gene_coordinates = dplyr::filter(hg38_gene_coordinates, gene_name %in% gene_symbol)
    }

    if(missing(gene_symbol) && !missing(ensembl_id)){
      gene_coordinates = dplyr::filter(hg38_gene_coordinates, ensembl_gene_id %in% gene_symbol)
    }
  }

  if(return_as == "bed"){
    region = paste0(gene_coordinates$chromosome, ":", gene_coordinates$start, "-", gene_coordinates$end)
  }else if(return_as == "df"){
    region = dplyr::select(gene_coordinates, chromosome, start, end, gene_name, hugo_symbol, ensembl_gene_id) %>%
      as.data.frame() %>%
      dplyr::arrange(chromosome, start)
  }

  if(!missing(gene_symbol)){
    message(paste0(nrow(region), " region(s) returned for ", length(gene_symbol), " gene(s)"))
  }

  if(!missing(ensembl_id)){
    message(paste0(nrow(region), " region(s) returned for ", length(ensembl_id), " gene(s)"))
  }

  return(region)
}

#' Return gennes residing in defined region(s)
#'
#' @param region Regions to intersect genes with, this should be a bed-like df (chromosome, start, end).
#' @param gene_format Parameter for specifying the format of returned genes, default is "hugo", other acceptable inputs are "ensembl".
#' @param genome_build Reference genome build.
#' @param chr_select Optional parameter to subset plot to specific chromosomes. Default value is chr1-22.
#'
#' @import data.table
#' @return
#' @export
#'
#' @examples
#' myc_region = gene_to_region(gene_symbol = "MYC", genome_build = "grch37", return_as = "df")
#' region = region_to_gene(region = myc_region, gene_format = "hugo", genome_build = "grch37")
#'
region_to_gene = function(region,
                          gene_format = "hugo",
                          genome_build = "grch37",
                          chr_select = paste0("chr", c(1:22))){

  if(genome_build == "grch37"){
    gene_list = grch37_gene_coordinates
  }

  if(genome_build == "hg38"){
    gene_list = hg38_gene_coordinates
  }

  #transform regions to data tables
  region_table = as.data.table(region)
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
    genes = select(inter_df, chromosome, start, end, hugo_symbol, region_start, region_end)
  }else if(gene_format == "ensembl"){
    genes = select(inter_df, chromosome, start, end, ensembl_gene_id, region_start, region_end)
  }

  #paste chr in chromosome column, if not there
  if(!str_detect(genes$chromosome, "chr")){
    genes = mutate(genes, chromosome = paste0("chr", chromosome))}

  genes = genes[genes$chromosome %in% chr_select, ]

  genes = as.data.frame(genes) %>%
      dplyr::arrange(chromosome, start)

  message(paste0(nrow(genes), " gene(s) returned for ", nrow(region), " region(s)"))

  return(genes)
}


#' Get a MAF that is just the variants unique to one of two flavours of variant calls available.
#'
#' @param these_sample_ids sample IDs to be included.
#' @param flavour1 First flavour of variant calls, to be returned as unique if not present in flavour2. Default is "clustered".
#' @param flavour2 Second flavour of variant calls.
#'
#' @return a list with MAFs that are only present in flavour1.
#' @export
#'
#' @examples
#' comp_flav = compare_mutation_flavour("clustered", "another_flavour")
#'
compare_mutation_flavour = function(these_sample_ids,
                                    flavour1 = "clustered",
                                    flavour2 = ""){

  these_dfs = list()
  for(this_sample_id in these_sample_ids){
    message(this_sample_id)
    maf1 = get_ssm_by_sample(this_sample_id, flavour = flavour1)
    maf2  = get_ssm_by_sample(this_sample_id, flavour = flavour2)
    maf1_only = intersect_maf(maf1, maf2)
    these_dfs[[this_sample_id]] = maf1_only
  }
  this_maf = rbindlist(these_dfs, use.names = TRUE)
  return(this_maf)
}


#' Perform set operations on two MAFs.
#'
#' @param maf1 First list of MAFs.
#' @param maf2 Second list of MAFs.
#' @param set_returned List of MAFs that doesn't share the same start positions as the other list of MAFs. Accepted commands are; "maf1_only" and "maf2_only", default is "maf1_only".
#'
#' @return Set of MAFs with start positions that don't match the start positions in the other supplied MAF file.
#' @export
#'
#' @examples
#' intersected_mafs_l1 = intersect_maf(maf_list1, maf_list2, "maf1_only")
#' intersected_mafs_l2 = intersect_maf(maf_list1, maf_list2, "maf2_only")
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


#' Tabulate mutation status for non-silent SSMs for a set of genes.
#'
#' @param gene_symbols List of gene symbols for which the mutation status will be tabulated. If not provided, lymphoma genes will be returned by default.
#' @param these_samples_metadata The matedata for samples of interest to be included in the returned matrix. Only the column "sample_id" is required. If not provided, the matrix is tabulated for all available samples as default.
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations. Default is TRUE.
#' @param augmented default: TRUE. Set to FALSE if you instead want the original MAF from each sample for multi-sample patients instead of the augmented MAF.
#' @param min_read_support Only returns variants with at least this many reads in t_alt_count (for cleaning up augmented MAFs).
#' @param maf_path If the status of coding SSM should be tabulated from a custom maf file, provide path to the maf in this argument. The default is set to NULL.
#' @param maf_data Either a maf loaded from disk or from the database using a get_ssm function.
#' @param include_hotspots Logical parameter indicating whether hotspots object should also be tabulated. Default is TRUE.
#' @param recurrence_min Integer value indicating minimal recurrence level.
#' @param review_hotspots Logical parameter indicating whether hotspots object should be reviewed to include functionally relevant mutations or rare lymphoma-related genes. Default is TRUE.
#' @param genes_of_interest List of genes for hotspot review. Currently only FOXO1, MYD88, and CREBBP are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF file. The default is hg19 genome build.
#' @param include_silent Logical parameter indicating whether to include siment mutations into coding mutations. Default is TRUE.
#'
#' @return Data frame.
#' @export
#'
#' @examples
#' coding_tabulated_df = get_coding_ssm_status(maf_data = "maf", gene_symbols=c("MYC","KMT2D"))
#' coding_tabulated_df = get_coding_ssm_status() #all lymphoma genes from bundled NHL gene list
#'
get_coding_ssm_status = function(gene_symbols,
                                 these_samples_metadata,
                                 from_flatfile = TRUE,
                                 augmented = TRUE,
                                 min_read_support = 3,
                                 maf_path = NULL,
                                 maf_data,
                                 include_hotspots = TRUE,
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
    }
  }
  return(all_tabulated)
}


#' INTERNAL FUNCTION called by prettyOncoplot, not meant for out-of-package usage.
#' Unsure what this function actually does, old helper function that is outdated?
#'
#' @param x
#'
#' @return Numeric value.
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


#' Count the number of mutations in a sliding window across a region for all samples. Unlikely to be used directly in most cases. See get_mutation_frequency_bin_matrix instead.
#'
#' @param this_region Genomic region in bed format.
#' @param chromosome Chromosome name in region.
#' @param start_pos Start coordinate of region.
#' @param end_pos End coordinate of region.
#' @param metadata Data frame containing sample ids and column with annotated data for the 2 groups of interest. All other columns are ignored. Currently, function exits if asked to compare more than 2 groups.
#' @param slide_by Slide size for sliding window, default is 100.
#' @param window_size Size of sliding window, default is 1000.
#' @param plot_type Set to TRUE for a plot of your bins. By default no plots are made.
#' @param sortByColumns Which of the metadata to sort on for the heatmap
#' @param return_format Return format of mutations. Accepted inputs are "long" and "long-simple". Default is "long-simple".
#' @param min_count_per_bin Minimum counts per bin, default is 3.
#' @param return_count Boolean statement to return count. Default is FALSE.
#' @param drop_unmutated This may not currently work properly. Default is FALSE. Not a very good parameter desciption?
#' @param classification_column Only used for plotting, default is "lymphgen"
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details). Default is FALSE.
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#' @return Count matrix.
#' @export
#' @import dplyr data.table ggplot2 cowplot
#'
#' @examples
#' chr11_mutation_freq = calc_mutation_frequency_sliding_windows(this_region = "chr11:69455000-69459900", metadata = meta_df, slide_by = 10, window_size = 10000)
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
                                                   mode = "slms-3",
                                                   ssh_session=NULL){

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
  region_ssm = GAMBLR::get_ssm_by_region(region = this_region, streamlined = FALSE, seq_type=seq_type, from_indexed_flatfile = from_indexed_flatfile, mode = mode,ssh_session=ssh_session) %>%
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


#' Write bedpe format data frame to a file that will work with IGV.
#'
#' @param sv_df data frame of bedpe formatted SV data.
#' @param filename File name (will be written to results/icgc_dart/misc/FILENAME). Default is "something.bedpe". Maybe not a very good default name?
#' @param add_chr_prefix Whether to force chr to be added to chromosome names. Default is TRUE.
#'
#' @return bedpe data frame that is compatible with IGN browser.
#' @export
#'
#' @examples
#' SVs_bedpe = sv_to_bedpe_file(sv_dataframe, "SVs.bedpe", TRUE)
#'
sv_to_bedpe_file = function(sv_df,
                            filename = "something.bedpe",
                            add_chr_prefix = TRUE){

  #add chr prefix if missing
  if(add_chr_prefix){
    if(!any(grepl("chr", region_sv$CHROM_A[1]))){
      sv_df = mutate(sv_df, CHROM_A = paste0("chr", CHROM_A)) %>%
        mutate(CHROM_B = paste0("chr", CHROM_B))
    }
  }
  bed_file = paste0("/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/misc/", filename)
  write.table(sv_df, file = bed_file, sep = "\t", quote = F, row.names = F, col.names = F)
}

#' INTERNAL FUNCTION called by calc_mutation_frequency_sliding_windows, not meant for out-of-package usage.
#' Parse a region string into chomosome, start and end.
#'
#' @param region A region string e.g. "chrX:12345-678910".
#'
#' @return A named list.
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


#' Convert mutation data to a shereable format.
#'
#' `sanitize_maf_data` returns oncomatrix of patient/gene data indicating only data needed to produce oncoplot.
#'
#' Write an oncomatrix from a MAF File for further plotting. This is meant to be run by individuals who have access to data sets to
#' "sanitize" a subset of data for subsequent use by them or others who don't have permission to access the raw data.
#' Example: User J has full permissions for ICGC data and has read permissions on a MAF file. User B needs to make some oncoplots
#' and/or perform some statistical analysis on the frequency and assortment of mutations in that data set but doesn't need all the details.
#' User J can run this function on a maf file and provide the path of the output to user B.
#'
#' @param mutation_maf_path Provide either the full path to a MAF file.
#' @param mutation_maf_data Otherwise provide a data frame of the MAF data.
#' @param output_oncomatrix Optionally provide the path for your sanitized output file (otherwise it writes to working directory).
#' @param genes_keep Specify which genes you want to remain in the output.
#' @param genes_drop Optionally specify which genes to drop (this doesn't mean all other genes will remain. Maftools decides that part).
#'
#' @return The full path to the oncomatrix file (a matrix with Variant_Classification or Multi_Hit indicating coding mutation status per patient).
#' @export
#'
#' @examples
#' lymph_genes = lymphoma_genes$Gene #note, this will be used by default if the user wants to be lazy
#' secure_maf = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/slms-3_vcf2maf_current/level_3/final_merged_grch37.CDS.maf"
#' safe_oncomatrix_path = sanitize_maf_data(mutation_maf_path = secure_maf, genes_keep = lymph_genes)
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
    genes_keep = lymphoma_genes$Gene
  }
  maf_o = maftools::read.maf(mutation_maf_data)

  maftools::oncoplot(maf_o, genes = genes_keep, writeMatrix = T, removeNonMutated = F)  #writes to working directory
  if(!missing(output_oncomatrix)){
    #rename it
    file.dplyr::rename("onco_matrix.txt", output_oncomatrix)
  }else{
    output_oncomatrix = paste0(getwd(), "/onco_matrix.tsv")
  }
  message(paste("your data is in:", output_oncomatrix))
  return(output_oncomatrix)
}


#' Annotate MAF-like data frome with a hot_spot column indicating recurrent mutations.
#'
#' @param mutation_maf A data frame in MAF format.
#' @param recurrence_min minimum number of recurrences for mutation to be included, default is 5.
#' @param analysis_base Base name for hot spot output directory.
#' @param p_thresh P value threshold, default is 0.05.
#'
#' @return The same data frame with one additional column "hot_spot".
#' @export
#'
#' @examples
#' hot_ssms = annotate_hotspots(all_ssm)
#' hot_maf = read.maf(hot_ssms)
#' oncoplot(hot_maf,genes=c("MEF2B","TP53","MYD88"),additionalFeature = c("hot_spot",TRUE))
#'
annotate_hotspots = function(mutation_maf,
                             recurrence_min = 5,
                             analysis_base = c("FL--DLBCL", "BL--DLBCL"),
                             p_thresh = 0.05){

  hotspot_info = list()
  for(abase in analysis_base){
    base_path = config::get("repo_base")

    clust_full_path = paste0(base_path, config::get("results_versioned")$oncodriveclustl$clusters)
    clust_full_path = glue::glue(clust_full_path)
    all_full_path = paste0(base_path, config::get("results_versioned")$oncodriveclustl$elements)
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


#' Annotate MAF-like data frome with a hot_spot column indicating recurrent mutations.
#'
#' @param annotated_maf A data frame in MAF format that has hotspots annotated using function annotate_hotspots().
#' @param genes_of_interest List of genes for hotspot review. Currently only FOXO1, MYD88, and CREBBP are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF file. The default is hg19 genome build.
#'
#' @return The same data frame with reviewed column "hot_spot".
#' @export
#' @import dplyr
#'
#' @examples
#' hot_ssms = review_hotspots(annotate_hotspots(get_coding_ssm()), genes_of_interest = c("CREBBP"))
review_hotspots = function(annotated_maf,
                          genes_of_interest = c("FOXO1", "MYD88", "CREBBP", "NOTCH1", "NOTCH2", "CD79B", "EZH2"),
                          genome_build = "hg19"){

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


#' Make a UCSC-ready custom track file from SV data.
#
#' @param sv_bedpe A bedpe formatted data frame of SVs.
#' @param output_file A bed file with UCSC custom header.
#' @param is_annotated Set to TRUE if input SV bedpe is annotated, default is TRUE.
#' @param sv_name SV name. Default is set to "all" = include all subtypes of SVs.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' all_sv = get_manta_sv()
#' sv_to_custom_track(all_sv, output_file = "GAMBL_sv_custom_track.bed", is_annotated = FALSE)
#'
sv_to_custom_track = function(sv_bedpe,
                              output_file,
                              is_annotated = TRUE,
                              sv_name = "all"){

  if(is_annotated){
  #reduce to a bed-like format
  sv_data1 = mutate(sv_bedpe, annotation = paste0(chrom1, ":", start1, "_", fusion)) %>%
    dplyr::select(chrom2, start2, end2, tumour_sample_id, annotation, fusion)

  sv_data2 = mutate(sv_bedpe, annotation = paste0(chrom2, ":", start2, "_", fusion)) %>%
    dplyr::select(chrom1, start1, end1, tumour_sample_id, annotation, fusion)

  print(head(sv_data1))
  print(head(sv_data2))
  colnames(sv_data1) = c("chrom", "start", "end", "sample_id", "annotation", "fusion")
  colnames(sv_data2) = c("chrom", "start", "end", "sample_id", "annotation", "fusion")
  sv_data = bind_rows(sv_data1, sv_data2)
  sv_data = mutate(sv_data, end = end + 10)
  }else{
    colnames(sv_data)[c(1,2,3)]=c("CHROM_A" ,  "START_A", "END_A" )
    colnames(sv_data)[c(4,5,6)]=c("CHROM_B" ,  "START_B", "END_B" )
    
    sv_data_1 = mutate(sv_bedpe, annotation = paste0( CHROM_B, ":", START_B)) %>%
      dplyr::select(CHROM_A, START_A, END_A, tumour_sample_id, annotation)
    sv_data_2 = mutate(sv_bedpe, annotation = paste0( CHROM_A, ":", START_A)) %>%
      dplyr::select(CHROM_B, START_B, END_B, tumour_sample_id, annotation)
    colnames(sv_data_1)=c("chrom", "start", "end", "sample_id", "annotation")
    colnames(sv_data_2)=c("chrom", "start", "end", "sample_id", "annotation")
    sv_data= bind_rows(sv_data_1,sv_data_2)
    
   # sv_data = dplyr::select(sv_data,chrom,start,end,sample_id,annotation)
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


#' Convert a maf-formatted data frame into a bed custom track file for UCSC.
#'
#' @param maf_data Either a maf loaded from disk or from the database using a get_ssm function.
#' @param output_file Name for your new bed file that can be uploaded as a custom track to UCSC.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' maf_to_custom_track(my_maf_data, "/home/rmorin/private/some_mutations.bed")
#'
maf_to_custom_track = function(maf_data,
                               these_samples_metadata,
                               output_file,
                               track_name="GAMBL mutations",
                               track_description = "mutations from GAMBL"){

  #reduce to a bed-like format
  maf_data = dplyr::select(maf_data, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)
  colnames(maf_data) = c("chrom", "start", "end", "sample_id")
  if(!any(grepl("chr", maf_data[,1]))){
    #add chr
    maf_data[,1] = unlist(lapply(maf_data[,1], function(x){paste0("chr", x)}))
  }
  lymphgen_cols = get_gambl_colours("lymphgen")
  colour_df = data.frame(lymphgen = names(lymphgen_cols), colour = lymphgen_cols)
  rgb_df = data.frame(t(col2rgb(lymphgen_cols))) %>%
    mutate(lymphgen = names(lymphgen_cols)) %>%
    unite(col = "rgb", red, green, blue, sep = ",")
  if(missing(these_samples_metadata)){
  meta = get_gambl_metadata() %>%
    dplyr::select(sample_id, lymphgen)
  }else{
    meta = these_samples_metadata %>%
      dplyr::select(sample_id, lymphgen)
  }
  samples_coloured = left_join(meta, rgb_df)
  maf_bed = maf_data %>%
    mutate(score = 0, strand = "+", start1 = start-1,start=start1, end1 = end)

  maf_coloured = left_join(maf_bed, samples_coloured, by = "sample_id") %>%
    dplyr::select(-lymphgen) %>%
    dplyr::filter(!is.na(rgb))
  header_ucsc = paste0('track name="',track_name,'" description="', track_description, '" visibility=2 itemRgb="On"\n')
  cat(header_ucsc,file = output_file)
  tabular = write.table(maf_coloured, file = output_file, quote = F, sep = "\t", row.names = F, col.names = F, append = TRUE)
}

test_glue = function(placeholder="INSERTED"){
  some_string = "this text has {placeholder}"
  print(some_string)
  ss=glue::glue(some_string)
  print(ss)
}

#' Bring together all derived sample-level results from many GAMBL pipelines.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param write_to_file Boolean statement that outputs tsv file if TRUE, default is FALSE.
#' @param join_with_full_metadata Join with all columns of meta data, default is FALSE.
#' @param these_samples_metadata Optional argument to use a user specified metadata df, overwrites get_gambl_metadata in join_with_full_metadata.
#' @param case_set Optional short name for a pre-defined set of cases.
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param from_cache Boolean variable for using cached results, default is TRUE.
#'
#' @return A table keyed on biopsy_id that contains a bunch of per-sample results from GAMBL
#' @export
#' @import tidyverse config
#'
#' @examples
#' everything_collated = collate_results(join_with_full_metadata = TRUE)
#'
collate_results = function(sample_table,
                           write_to_file = FALSE,
                           join_with_full_metadata = FALSE,
                           these_samples_metadata,
                           case_set,
                           sbs_manipulation = "",
                           seq_type_filter = "genome",
                           from_cache = TRUE,
                           ssh_session){

  # important: if you are collating results from anything but WGS (e.g RNA-seq libraries) be sure to use biopsy ID as the key in your join
  # the sample_id should probably not even be in this file if we want this to be biopsy-centric
  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter = seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  if(write_to_file){
    from_cache = FALSE #override default automatically for nonsense combination of options
  }
  if(from_cache){
    output_file = config::get("table_flatfiles")$derived
    output_base = config::get("repo_base")
    output_file = paste0(output_base, output_file)
    output_file = glue::glue(output_file)

    #check for missingness
    if(!file.exists(output_file)){
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    sample_table = read_tsv(output_file) %>% dplyr::filter(sample_id %in% sample_table$sample_id)
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
    output_file = config::get("table_flatfiles")$derived
    output_file = glue::glue(output_file)
    output_base = config::get("repo_base")
    output_file = paste0(output_base, output_file)
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


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#' Extract derived results stored in the database (these are usually slower to derive on the fly).
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations. Default is FALSE, TRUE is not yet implemented.
#'
#' @return Data frame with one row per sample. Contains the contents of the derived_data table in the database.
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' gambl_results_derived = collate_derived_results(samples_df)
#'
collate_derived_results = function(sample_table,
                                   seq_type_filter="genome",
                                   from_flatfile = FALSE){

  if(from_flatfile){
    message("not implemented YET")
  }else{
    database_name = config::get("database_name")
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


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#' Collate a few CSR annotations, including MiXCR.
#'
#' @param sample_table A data frame with sample_id as the first column.
#'
#' @return The sample table with additional columns.
#' @import tidyverse
#'
#' @examples
#' gambl_results_derived = collate_csr_results(gambl_results_derived)
#'
collate_csr_results = function(sample_table,seq_type_filter = "genome"){
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


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#' Compute some summary statistics based on SSM calls.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations, default is TRUE.
#' @param include_silent Logical parameter indicating whether to include siment mutations into coding mutations. Default is FALSE.
#'
#' @return Sample table.
#'
#' @examples
#' ssm_results = colalte_ssm_results(samples, TRUE, TRUE)
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
    base_path = config::get("project_base")
    #test if we have permissions for the full gambl + icgc merge
    maf_partial_path = config::get("results_flatfiles")$ssm$template$merged$deblacklisted

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


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#' Collate all SV calls from the genome data and summarize for main oncogenes of interest per sample.
#'
#' @param sample_table A data frame with sample_id as the first column.
#'
#' @return The sample table with additional columns.
#' @import tidyverse
#'
#' @examples
#' gambl_results_derived = collate_curated_sv_results(gambl_results_derived)
#'
collate_curated_sv_results = function(sample_table,seq_type_filter="genome"){

  path_to_files = config::get("derived_and_curated")
  project_base = config::get("project_base")
  manual_files = dir(paste0(project_base, path_to_files), pattern = ".tsv")
  for(f in manual_files){
    full = paste0(project_base, path_to_files, f)
    this_data = suppressMessages(read_tsv(full, comment = "#"))
    #TO DO: fix this so it will join on biopsy_id or sample_id depending on which one is present, Done?
    sample_table = left_join(sample_table, this_data)
  }
  return(sample_table)
}

#' Get wildcards for a sample_id/seq_type combination.
#'
#' @param this_sample_id
#' @param seq_type
#'
#' @return
#' @export
#'
#' @examples
get_sample_wildcards = function(this_sample_id,seq_type){
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

#' Annotate mutations with their copy number information.
#'
#' @param this_sample Sample ID of the sample you want to annotate.
#' @param coding_only Optional. set to TRUE to rescrict to only coding variants.
#' @param from_flatfile Optional. Instead of the database, load the data from a local MAF and seg file.
#' @param use_augmented_maf Boolean statement if to use augmented maf, default is FALSE.
#' @param maf_file Path to maf file.
#' @param seq_file path to seq file.
#' @param seg_file_source Specify what copy number calling program the input seg file is from, as it handles ichorCNA differently than WisecondorX, Battenberg, etc.
#' @param assume_diploid Optional. If no local seg file is provided, instead of defaulting to a GAMBL sample, this parameter annotates every mutation as copy neutral.
#' @param genes Genes of interest.
#' @param include_silent Logical parameter indicating whether to include siment mutations into coding mutations. Default is FALSE
#'
#' @return A list containing a data frame (MAF-like format) with two extra columns:
#' log.ratio is the log ratio from the seg file (NA when no overlap was found)
#' as well as the segmented copy number data with the same copy number information
#' CN is the rounded absolute copy number estimate of the region based on log.ratio (NA when no overlap was found)
#' @export
#' @import tidyverse data.table RMariaDB DBI dbplyr
#'
#' @examples
#' cn_list = assign_cn_to_ssm(this_sample = "HTMCP-01-06-00422-01A-01D", coding_only = TRUE)
#'
assign_cn_to_ssm = function(this_sample,
                            coding_only = FALSE,
                            from_flatfile = TRUE,
                            use_augmented_maf = TRUE,
                            tool_name = "battenberg",
                            maf_file,
                            seg_file,
                            seg_file_source = "battenberg",
                            assume_diploid = FALSE,
                            genes,
                            include_silent = FALSE,
                            ssh_session,
                            seq_type="genome",
                            projection="grch37"){

  database_name = config::get("database_name")
  project_base = config::get("project_base")
  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  if(!missing(maf_file)){
    maf_sample = fread_maf(maf_file) %>%
      dplyr::mutate(Chromosome = gsub("chr", "", Chromosome))
  }else if(from_flatfile){
    #get the genome_build and other wildcards for this sample
    wildcards = get_sample_wildcards(this_sample,seq_type)
    genome_build = wildcards$genome_build
    unix_group = wildcards$unix_group
    seq_type = wildcards$seq_type
    tumour_sample_id = wildcards$tumour_sample_id
    normal_sample_id = wildcards$normal_sample_id
    pairing_status = wildcards$pairing_status
    maf_sample = get_ssm_by_sample(this_sample_id = this_sample, augmented = use_augmented_maf, ssh_session = ssh_session)

  }else{
    #get all the segments for a sample and filter the small ones then assign CN value from the segment to all SSMs in that region
    con = dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    maf_table = config::get("results_tables")$ssm
    maf_sample = dplyr::tbl(con, maf_table) %>%
      dplyr::filter(Tumor_Sample_Barcode == this_sample) %>%
      as.data.frame()
  }
  if(coding_only){
    maf_sample = dplyr::filter(maf_sample, Variant_Classification %in% coding_class)
  }

  if(!missing(genes)){
    maf_sample = dplyr::filter(maf_sample, Hugo_Symbol %in% genes)
  }

  if(!missing(seg_file)){
    seg_sample = read_tsv(seg_file) %>%
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
    project_base = config::get("project_base",config="default")
    local_project_base = config::get("project_base")

    results_path_template = config::get("results_flatfiles")$cnv$battenberg
    results_path = paste0(project_base, results_path_template)
    local_results_path = paste0(local_project_base, results_path_template)

    ## NEED TO FIX THIS TO contain tumour/normal ID from metadata and pairing status
    battenberg_file = glue::glue(results_path)
    local_battenberg_file = glue::glue(local_results_path)


    message(paste("looking for flatfile:", battenberg_file))
    if(!missing(ssh_session)){
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
      message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
      message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
    }

    seg_sample = read_tsv(battenberg_file) %>%
      as.data.table() %>%
      dplyr::mutate(size = end - start) %>%
      dplyr::filter(size > 100) %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
      dplyr::rename(Chromosome = chrom, Start_Position = start, End_Position = end)

    data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    a = data.table::as.data.table(maf_sample)
  }else{
    seg_sample = get_sample_cn_segments(this_sample_id = this_sample) %>%
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


#' Estimate purity.
#'
#' @param in_maf Path to a local maf file.
#' @param in_seg Path to a local corresponding seg file for the same sample ID as the input maf.
#' @param sample_id Specify the sample_id or any other string you want embedded in the file name.
#' @param seg_file_source Specify what copy number calling program the input seg file is from, as it handles ichorCNA differently than WisecondorX, battenberg, etc.
#' @param show_plots Optional. Show two faceted plots that display the VAF and purity distributions for each copy number state in the sample.
#' @param assume_diploid Optional. If no local seg file is provided, instead of defaulting to a GAMBL sample, this parameter annotates every mutation as copy neutral.
#' @param coding_only Optional. set to TRUE to rescrict to only coding variants.
#' @param genes Genes of interest.
#'
#' @return A list containing a data frame (MAF-like format) with the segmented absolute copy number data and three extra columns:
#' VAF is the variant allele frequency calculated from the t_ref_count and t_alt_count
#' Ploidy is the number of copies of an allele in the tumour cell
#' Final_purity is the finalized purity estimation per mutation after considering different copy number states and LOH events
#' @export
#' @import tidyverse
#'
#' @examples
#' tumour_purity = estimate_purity(in_maf = "path/to/file.maf", in_seg = "path/to/file.seg", seg_file_source = "ichorCNA", show_plots = TRUE)
#' tumour_purity = estimate_purity(in_maf="path/to/file.maf", assume_diploid = TRUE)
#'
estimate_purity = function(in_maf,
                           in_seg,
                           sample_id,
                           seg_file_source = "battenberg",
                           show_plots = FALSE,
                           assume_diploid = FALSE,
                           coding_only = FALSE,
                           genes,
                           ssh_session){

  # Merge the CN info to the corresponding MAF file, uses GAMBLR function
  if(missing(in_maf) & missing(in_seg)){
    CN_new = assign_cn_to_ssm(this_sample = sample_id, coding_only = coding_only, assume_diploid = assume_diploid, genes = genes,seg_file_source = seg_file_source,ssh_session=ssh_session)$maf
  }else if(!missing(in_seg)){
    CN_new = assign_cn_to_ssm(this_sample = sample_id, maf_file = in_maf, seg_file = in_seg, seg_file_source = seg_file_source, coding_only = coding_only, genes = genes,ssh_session=ssh_session)$maf
  }else{
    # If no seg file was provided and assume_diploid paramtere is set to true,
    if(assume_diploid){
      CN_new = assign_cn_to_ssm(this_sample = sample_id, maf_file = in_maf, assume_diploid = TRUE, coding_only = coding_only, genes = genes,ssh_session=ssh_session)$maf
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


#' INTERNAL FUNCTION called by referesh_metadata_tables, not meant for out-of-package usage.
#' Refresh the contents of a database table.
#'
#' @param table_name Name of table to refresh.
#' @param connection Database connection object.
#' @param file_path Path to the table contents to populate.
#'
#' @return Table.
#'
#' @examples
#' refresh_full_table(table_x, con,file_x)
#'
refresh_full_table = function(table_name,
                              connection,
                              file_path){

  table_data = read_tsv(file_path)
  dbWriteTable(con, table_name, table_data, overwrite = TRUE)
  print(paste("POPULATING table:", table_name, "USING path:", file_path))
}


#' INTERNAL FUNCTION, not meant for out-of-package usage.
#' Refresh the contents of meta data table.
#'
#' @return Table.
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


#' Function that performes sanity checks on meta data.
#'
#' @return Table.
#' @export
#'
#' @examples
#' sane_meta_data = sanity_check_metadata()
#'
sanity_check_metadata = function(){

  cfg = config::get("tables")
  database_name = config::get("database_name")
  metadata_tables = tibble(key = names(cfg), table = cfg) %>%
    unnest_auto("table")

  cfg = config::get("table_flatfiles")
  metadata_files = tibble(key = names(cfg), file = cfg) %>%
    unnest_auto("file")

  all_metadata_info = left_join(metadata_tables, metadata_files)
  base_path = config::get("repo_base")
  all_metadata_info = all_metadata_info %>%
    mutate(file = paste0(base_path, file))

  all_metadata_df = all_metadata_info %>%
    column_to_rownames(var = "key")
  #all samples with different seq_type and protocol must have a unique sample_id
  sample_df = read_tsv(all_metadata_df["samples", "file"])
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


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param somalier_output Somalier ancestery.tsv
#'
#' @return Table.
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
  somalier_all = read_tsv(somalier_output)
  somalier_all = mutate(somalier_all, sample_id = str_remove(`#sample_id`, pattern = ":.+")) %>%
    dplyr::select(-`#sample_id`, -given_ancestry)
  somalier_all = dplyr::select(somalier_all, sample_id, predicted_ancestry, PC1, PC2, PC3, PC4, PC5)
  sample_table = left_join(sample_table, somalier_all)
  return(sample_table)
}


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param file_path Path to extra emtadata.
#'
#' @return Table.
#'
#' @examples
#' table = collate_extra_metadata(sample_table = "my_sample_table.txt")
#'
collate_extra_metadata = function(sample_table,
                                  file_path){

  file_path = "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/experiments/2021-04-21-Trios-MiXCR/trios_relapse_timing.tsv"
  extra_df = read_tsv(file_path)
  sample_table = left_join(sample_table, extra_df, by = c("sample_id" = "biopsy_id"))
}


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#' Bring in the results from mutational signature analysis.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param file_path Optional path to SBS file.
#' @param scale_vals Parameter not used?
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#'
#' @return A data frame with new columns added.
#' @import tidyverse
#'
#' @examples
#' collated = collate_sbs_results(sample_table=sample_table,sbs_manipulation=sbs_manipulation)
#'
collate_sbs_results = function(sample_table,
                               seq_type_filter="genome",
                               file_path,
                               scale_vals = FALSE,
                               sbs_manipulation = ""){
  if(seq_type_filter!="genome"){
    message("skipping sbs for seq_type")
    return(sample_table)
  }
  if(missing(file_path)){
    base = config::get("project_base")

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


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#' Determine which cases have NFKBIZ UTR mutations.
#'
#' @param sample_table A data frame with sample_id as the first column.
#'
#' @return Samples table.
#' @import tidyverse
#'
#' @examples
#' sample_table = collate_nfkbiz_results(sample_table=sample_table)
#'
collate_nfkbiz_results = function(sample_table,seq_type_filter="genome",ssh_session=NULL){
  #TO DO: Update to work with hg38 projection
  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter=seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  this_region = "chr3:101578214-101578365"
  nfkbiz_ssm = get_ssm_by_region(region = this_region,seq_type = seq_type_filter,ssh_session=ssh_session) %>%
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


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#' Determine the hypermutation status of a few genes.
#'
#' @param sample_table A data frame with sample_id as the first column.
#'
#' @return Samples table.
#' @import tidyverse
#'
#' @examples
#' sample_table = collate_ashm_results(sample_table=sample_table)
#'
collate_ashm_results = function(sample_table,
                                seq_type_filter="genome",
                                ssh_session=NULL){

  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter=seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  #just annotate BCL2, MYC and CCND1 hypermutation
  regions_df = data.frame(name = c("CCND1","BCL2","MYC"),
  region = c("chr11:69455000-69459900", "chr18:60983000-60989000", "chr8:128747615-128751834"))
  region_mafs = lapply(regions_df$region, function(x){get_ssm_by_region(region = x, streamlined = FALSE,seq_type=seq_type_filter,ssh_session=ssh_session)})
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


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#' Determine and summarize which cases have specific oncogene SVs.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param tool Name of tool (optional, default is manta).
#' @param oncogenes Which oncogenes to collate SVs from.
#'
#' @return Data frame with additional columns ({tool}_{oncogene} and {tool}_{oncogene}_{partner}).
#' @import tidyverse
#'
#' @examples
#' results = collate_samples_sv_results(samples, "manta", c("MYC", "BCL2"))
#'
collate_sv_results = function(sample_table,
                              tool = "manta",
                              seq_type_filter="genome",
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


#' Get some colour schemes for annotating figures.
#'
#' @param classification (optionally request only colours for pathology, lymphgen, mutation or copy_number).
#' @param alpha alpha of plotted colours.
#'
#' @return A named vector of colour codes for lymphgen classes and pathology
#' @export
#' @import tidyverse
#'
#' @examples
#' lymphgen_cols = get_gambl_colours("lymphgen")
#' # be sure to install ggsci from https://github.com/morinlab/ggsci
#' # install_github("morinlab/ggsci")
#'
get_gambl_colours = function(classification = "all",
                             alpha = 1){

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
                          "DLBCL-A" = "#721F0F",
                          "IC-BL" = "#45425A",
                          "DGG-BL" = "#E90C8BFF",
                          "DLBCL-B" = "#FB9A99",
                          "DLBCL-C" = "#C41230")

  all_colours[["FL"]] = c(dFL = "#99C1B9", cFL = "#D16666", DLBCL = "#479450")
  
  all_colours[["lymphgenerator"]] = c("MP3"="#5B8565",
                                      "EGB" = "#98622A",
                                      "ETB"="#813F3D",
                                      "aSCI"="#D66B1F",
                                      "aSEL"="#C41230",
                                      "MCaP"="#5F8CFF",
                                      "BNZ"="#8870B6"
                                      )
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
        "3'UTR" = unname(blood_cols["Yellow"]))

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
      "CNS" = "#E2EF60",
      "THRLBCL" = "#A5F2B3",
      "MM"="#CC9A42",
      "SCBC"="#8c9c90",
      "UNSPECIFIED"="#cfba7c",
      "OTHER"="#cfba7c",
      "MZL"="#065A7F",
      "SMZL"="#065A7F"
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

  #print(all_colours)
  for(colslot in names(all_colours)){
    raw_cols = all_colours[[colslot]]
    raw_cols_rgb = col2rgb(raw_cols)
    alpha_cols = rgb(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ], alpha = alpha * 255L, names = names(raw_cols), maxColorValue = 255L)
    names(alpha_cols) = names(raw_cols)
    all_colours[[colslot]] = alpha_cols
  }
  for(this_group in names(all_colours)){
    everything = c(everything, all_colours[[this_group]])
  }
  #return matching value from lowercase version of the argument if it exists
  lc_class = stringr::str_to_lower(classification)
  if(classification %in% names(all_colours)){
    return(all_colours[[classification]])
  }else if(lc_class %in% names(all_colours)){
    return(all_colours[[lc_class]])
  }else{
    return(everything)
  }
}


#' Get full paths for bam files for a sample or patient.
#'
#' @param sample Either provide sample_id or patient_id.
#' @param patient Either provide sample_id or patient_id.
#'
#' @return A list that contains the genome_build and an igv-friendly build (igv_build), a list of bam file paths for tumour, normal and mrna data.
#' @export
#' @import tidyverse
#'
#' @examples
#'
#' this_sv = filter(annotate_sv(get_manta_sv()),partner=="HIST1H2BK")
#' #arbitrarily grab one SV
#' bam_details = get_bams(sample=this_sv[1,"tumour_sample_id"])
#'
get_bams = function(sample,
                    patient
                    ){

  meta = get_gambl_metadata(tissue_status_filter = c("tumour", "normal"),seq_type_filter="genome")
  meta_mrna = get_gambl_metadata(seq_type_filter = "mrna")
  #get all samples for this patient
  if(missing(patient)){
    patient = meta %>%
      dplyr::filter(sample_id == sample) %>%
      dplyr::pull(patient_id)
  }
  meta_patient = meta %>%
    dplyr::filter(patient_id == patient)

  meta_mrna_patient = meta_mrna %>%
    dplyr::filter(patient_id == patient)

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

  bam_details$pairing_status = dplyr::filter(meta_patient, tissue_status == "tumour") %>%
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


#' Load bams and generate an IGV screenshot for one or more regions.
#'
#' @param bams Character vector containing the full path to one or more bam files.
#' @param genome_build String specifying the genome build for the bam files provided.
#' @param region Optionally specify the region as a single string (e.g. "chr1:1234-1235").
#' @param padding Optionally specify a positive value to broaden the region around the specified position. Default is 200.
#' @param chrom Optionally specify the region by specifying the chromosome, start and end (see below).
#' @param start Optionally specify the region by specifying the start.
#' @param end Optionally specify the region by specifying the end.
#' @param sample_id Specify the sample_id or any other string you want embedded in the file name.
#' @param out_path Specify the output directory where the snapshot will be written.
#' @param igv_port Specify the port IGV is listening on.
#'
#' @return Path to file (.png)
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
#'
make_igv_snapshot = function(bams,
                             genome_build,
                             region,
                             padding = 200,
                             chrom,
                             start,
                             end,
                             sample_id,
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
  filename = paste(sample_id, region, "snapshot.png", sep = "_")
  IGVsnapshot(sock, fname = filename, dirname = out_path)
  return(paste0(out_path, filename))
}


#' Using GISTIC2.0 outputs, perform Fisher's exact test to compare CNV frequencies between 2 groups.
#'
#' @param gistic_lesions Path to the GISTIC2.0 all_lesions output file.
#' @param metadata Data frame containing sample ids and column with annotated data for the 2 groups of interest. All other columns are ignored. Currently, function exits if asked to compare more than 2 groups.
#' @param comparison Specify column annotating groups of interest.
#' @param fdr.method FDR method to adjust p values. Uses stats::p.adjust function, and therefore accepts its method for FDR ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). By default, this function uses "fdr".
#' @param fdr.cutoff Specify FDR significance cut-off. By default, this function uses 0.1.
#' @param text_size Size of the text on the forest plot of differentially enriched CNV. Default text-size is 7.
#' @param blacklisted_regions Optionally, specify any descriptors (value from column `Descriptor` of GISTIC2.0 all_lesions output file) to filter out before amy comparisons are done. It is possible to specify list of multiple descriptors, for example, c("3p12.3", "12p13.2"). Default is NULL.
#'
#' @return list of 3 objects:
#' DISTINCT - set of CNV identified as significantly different between 2 groups;
#' CNV.EVENTS - conveniently reformatted set of events in both groups that can be used for downstream analyses; and
#' GRAPH - forest plot visualization of distinct CNV between 2 groups.
#' Output can be accessed by index (e.g.[3]), or by name (e.g.["GRAPH"])
#' @export
#' @import tidyverse stats metaviz
#'
#' @examples
#' # basic usage
#' FtestCNV(gistic_lesions = "path_to_GISTIC2.0_output/all_lesions.conf_{confidence_level}.txt", metadata = derived_data, comparison = "pathology")
#' # advanced usage
#' FtestCNV(gistic_lesions = "path_to_GISTIC2.0_output/all_lesions.conf_{confidence_level}.txt", metadata = derived_data,
#' comparison = "pathology", fdr.method = "bonferroni", fdr.cutoff = 0.05, blacklisted_regions = c("3p12.3", "12p13.2"))
#'
FtestCNV = function(gistic_lesions,
                     metadata,
                     comparison,
                     fdr.method = "fdr",
                     fdr.cutoff = 0.1,
                     text_size = 7,
                     blacklisted_regions = NULL){

  # get groups of interest for comparison
  GROUPS.TO.COMPARE = unique(metadata[,comparison])

  # quick check that only 2 groups are provided forr comparison
  if(length(GROUPS.TO.COMPARE) > 2){
    message("The current implementation of function only accepts 2 groups for comparison. You provided 3 groups.")
    message("Please modify metadata accordingly to compare only 2 groups.")
    message("Groups you provided are: ", paste(c(GROUPS.TO.COMPARE), collapse = ","))
    return(NULL)
  }
  # read lesions from gistic utput to collect event/case
  lesions = readr::read_tsv(gistic_lesions, col_names = TRUE) %>%
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
  GROUP1_vs_GROUP2$FDR = stats::p.adjust(GROUP1_vs_GROUP2$pVal, method = fdr.method)

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


#' If some samples are missing from the matrix, add them with filled in 0 as value and normalize their ordering for consistency.
#'
#' @param incoming_matrix A matrix or data frame that should be filled. Required parameter.
#' @param list_of_samples Vector specifying all desired samples to be present in the resulting matrix. Required parameter.
#' @param fill_in_values Value that will be used to fill in the matrix.
#' @param normalize_order Logical parameter specifying whether sample order should be according to the supplied list. Default is TRUE.
#' @param samples_in_rows Logical argument indicating whether samples are in rows or columns. Default assumes samples are in rows and columns are features.
#'
#' @return A data frame with maintained orientation (rows and columns) where samples from the supplied list are present and reordered according to the specified order.
#' @export
#'
#' @examples
#' partial_matrix = get_coding_ssm_status(these_samples_metadata = (get_gambl_metadata(case_set = "BL--DLBCL") %>% filter(pairing_status == "unmatched")), include_hotspots = FALSE)
#' complete_matrix = complete_missing_from_matrix(partial_matrix, get_gambl_metadata() %>% pull(sample_id))
#'
complete_missing_from_matrix = function(incoming_matrix,
                                        list_of_samples,
                                        fill_in_values = 0,
                                        normalize_order = TRUE,
                                        samples_in_rows = TRUE){

  # check for required arguments
  if (missing(incoming_matrix)){
      stop("Please provide initial matrix to fill.")
  }

  if (missing(list_of_samples)){
      stop("Please provide list of samples to complete the matrix and normalize order.")
  }

  # is samples are in columns, transpose the matrix so code below is generalizable
  if(!samples_in_rows){
    incoming_matrix = as.data.frame(incoming_matrix) %>%
      t()
  }

  matrix_with_all_samples = rbind(incoming_matrix, matrix(fill_in_values:fill_in_values, # populate matrix with all 0
                                                          length(setdiff(list_of_samples, rownames(incoming_matrix))), # how many rows
                                                          ncol(incoming_matrix), # how many columns
                                                          dimnames = list(setdiff(list_of_samples, rownames(incoming_matrix)), # name rows with sample IDs
                                                          colnames(incoming_matrix))) %>% # name columns with gene names
                                                            as.data.frame(.))

  # this is very helpful in clustering
  if(normalize_order){
    matrix_with_all_samples = matrix_with_all_samples[ order(match(rownames(matrix_with_all_samples), list_of_samples)),]
  }

  # transpose matrix back to the initial format supplied by user (samples in columns)
  if(!samples_in_rows){
    matrix_with_all_samples = as.data.frame(matrix_with_all_samples) %>%
      t()
  }
  return(matrix_with_all_samples)
}


#' Subset maf file to only features that would be available in the WEX data.
#'
#' @param maf Incoming maf object. Can be maf-like data frame or maftools maf object. Required parameter. Minimum columns that should be present are Chromosome, Start_Position, and End_Position.
#' @param custom_bed Optional argument specifying a path to custom bed file for covered regions. Must be bed-like and contain chrom, start, and end position information in first 3 columns. Other columns are disregarded if provided.
#' @param genome_build String indicating genome build of the maf file. Default is grch37, but can accept modifications of both grch37- and hg38-based duilds.
#' @param padding Numeric value that will be used to pad probes in WEX data from both ends. Default is 100. After padding, overlapping features are squished together.
#' @param chr_prefixed Is the data chr-prefixed or not? Default is FALSE.
#'
#' @return Aa data frame of maf-like object with same columns as in input, but where rows are only kept for features that would be present as if the sample is WEX.
#' @export
#' @import tidyverse data.table
#'
#' @examples
#' myc_ashm_maf = get_ssm_by_region(region="8:128748352-128749427") # get all ssm in the MYC aSHM region
#' GAMBLR::genome_to_exome(maf = myc_ashm_maf) # get mutations with 100 bp padding (default)
#' GAMBLR::genome_to_exome(maf = myc_ashm_maf, padding = 0) # get mutations covered in WEX with no padding
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


#' Consolidate a column of LymphGen data in the original Subtype.Prediction output format to the GAMBLR tidy format
#'
#' @param df Input data frame.
#' @param lymphgen_column_in The name of the column with lymphgen data to be processed.
#' @param lymphgen_column_out The name of the column to write the tidied results (optional).
#' @param relevel If TRUE, will return the output column as a factor with plot-friendly levels.
#'
#' @return A data frame with a tidied lymphGen column
#' @export
#' @import tidyverse
#'
#' @examples
#' metadata <- get_gambl_metadata()
#' GAMBLR::tidy_lymphgen(metadata,  "lymphgen_with_cnv", "lymphgen_with_cnv_tidy", relevel = TRUE)
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
    str_detect(.data[[lymphgen_column_in]],"ST2")~"ST2-COMP",
    str_detect(.data[[lymphgen_column_in]],"BN2")~"BN2-COMP"
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

#' Supplement the "lymphgen" column of the metadata with classification for additional samples.
#' Expects at least to have columns "patient_id" to bind on, and "lymphgen" to supplement the data on.
#'
#'
#' @param sample_table Input data frame with metadata.
#' @param derived_data_path Optional argument specifying the path to a folder with files following the pattern *lymphgen.txt
#'
#' @return A data frame with a supplemented lymphGen column
#' @export
#' @import tidyverse
#'
#' @examples
#' metadata <- get_gambl_metadata()
#' GAMBLR::collate_lymphgen(metadata)
#'
collate_lymphgen = function(sample_table,
                            derived_data_path = "",
                            verbose = TRUE) {
  if (derived_data_path == "") {
    path_to_files = config::get("derived_and_curated")
    project_base = config::get("project_base")
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


#' INTERNAL FUNCTION called by collate_results, not meant for out-of-package usage.
#' Collate qc metrics.
#'
#' @param sample_table df with sample ids in the first column.
#' @param seq_type_filter default is genome, capture is also available for unix_group icgc_dart.
#'
#' @return The sample table with additional columns.
#' @import tidyverse
#'
#' @examples
#' sample_table <- data.frame (sample_id  = c("09-33003T", "05-18426T", "BLGSP-71-06-00174-01A-01D", "BLGSP-71-22-00347-01A-12E", "SP59282"))
#' qc_metrics = collate_qc_results(sample_table = sample_table, seq_type_filter = "genome")
#'
collate_qc_results = function(sample_table,
                              seq_type_filter = "genome"){

  if(! seq_type_filter %in% c("genome", "capture")){
    stop("Please provide a valid seq_type (\"genome\" or \"capture\").")
  }

  #get paths
  base = config::get("project_base") #base
  qc_template = config::get("qc_met") #qc metric template
  qc_path = glue::glue(qc_template) #glue wildcards
  icgc_qc_path_full = paste0(base, qc_path) #combine for icgc_dart qc metric path
  gambl_qc = gsub("icgc_dart", "gambl", qc_path) #use gsub to substitute icgc_dart in qc path with gambl
  gambl_qc_path_full = paste0(base, gambl_qc) #combine for gambl qc metric path

  #read in icgc qc data, rename sample id column and filter on samples in sample id in sample_table
  icgc_qc = suppressMessages(read_tsv(icgc_qc_path_full)) %>%
    rename(sample_id = UID) %>%
    dplyr::filter(sample_id %in% dplyr::pull(sample_table, sample_id))

  #read in gambl qc data (if seq_type_filter set to "genome"), rename sample id column and filter on samples in sample id in sample_table
  if(seq_type_filter == "genome"){
    gambl_qc = suppressMessages(read_tsv(gambl_qc_path_full)) %>%
      rename(sample_id = UID) %>%
      dplyr::filter(sample_id %in% dplyr::pull(sample_table, sample_id))

    #join gambl and icgc QC data
    full_qc = rbind(gambl_qc, icgc_qc)

    return(full_qc)

    }else{
      message("Currently, seq_type_filter = \"capture\" is only available for unix_group \"icgc_dart\". Only QC metrics for icgc_dart will be returned.")
      #TO DO: Remove this once capture metrics have been generated for gambl samples.
      return(icgc_qc)
      }
}


#' Will prepare the data frame of binary matrix to be used as NMF input. This means that for the features with SSM and CNV,
#' they will be squished together as one feature named GeneName-MUTorAMP or GeneName-MUTorLOSS, so the CNV features in the input data frame are expected
#' to be named GeneName_AMP or GeneName_LOSS. Next, for the genes with hotspot mutations labelled in the input data as
#' GeneNameHOTSPOT, the feature for hotspot mutation will be given preference and SSM with/without CNV will be set to 0 for that sample.
#' The naming scheme of the features as in this description is important, because the function uses regex to searh for these patters as specified.
#' Finally, if any features are provided to be dropped explicitly, they will be removed, and then the features not meeting the specified minimal
#' frequency will be removed, as well as any samples with 0 features.
#' Consistent with NMF input, in the input data frame each row is a feature, and each column is a sample. The input is expected to be numeric 1/0 with row and column names.
#'
#'
#' @param incoming_data Input data frame or matrix to prepare for NMF.
#' @param blacklisted_cnv_regex Regular expression to match in feature names when considering SSM/CNV overlap.
#' @param drop_these_features Optional argument with features to drop from resulting matrix.
#' @param min_feature_percent Minimum frequency for the feature to be returned in the resulting matrix. By default, features present in less than 0.5% of samples will be discarded.
#'
#' @return A matrix compatible with NMF input.
#' @export
#' @import tidyverse
#'
#' @examples
#' data = system.file("extdata", "sample_matrix.tsv", package = "GAMBLR") %>% read_tsv() %>% column_to_rownames("Feature")
#' NMF_input = massage_matrix_for_clustering(data)
#'

massage_matrix_for_clustering = function(incoming_data,
                                         blacklisted_cnv_regex="3UTR|SV|HOTSPOT|TP53BP1|intronic",
                                         drop_these_features,
                                         min_feature_percent = 0.005){

  # if there is a CNV and mutation at the same gene, squish these features together
  message("Searching for overlapping CNV and mutation features to squish together ...")
  feat_with_cnv_data = rownames(incoming_data)[grepl("AMP|LOSS", rownames(incoming_data))]
  output_data = incoming_data

  for (g in feat_with_cnv_data){
    this_feature = unlist(strsplit(g, split='_', fixed=TRUE))
    red_features <- rownames(output_data)[grepl(this_feature[1], rownames(output_data))]
    red_features <- red_features[!grepl(blacklisted_cnv_regex, red_features)] # these features to be kept separately
    if(length(red_features)>1){
      message(paste0("Found redundant features for gene ", red_features[1], ", processing ..."))
      output_data[,output_data[red_features[2],]>0][red_features,][red_features[1],] = 1
      rownames(output_data)[rownames(output_data)==red_features[1]] = paste0(this_feature[1],
                                                               "-MUTor",
                                                               this_feature[2])
      output_data = output_data[!rownames(output_data) %in% red_features[2],]

    }
  }

  message("Success")

  # if there is a hotspot and SSM for same gene, give priority to hotspot
  message("Searching for overlapping HOTSPOT and mutation features to squish together ...")
  feat_with_hotspot_data = rownames(output_data)[grepl("HOTSPOT", rownames(output_data))]
  for (hot in feat_with_hotspot_data){
    this_gene=gsub("HOTSPOT","", hot)
    # this gene may also have CNV data already squished
    maybe_cnv = grepl("MUTor",
                      rownames(output_data[grepl(this_gene,
                                          rownames(output_data)),]))
    if("TRUE" %in% maybe_cnv){ # if it has the cnv data then use the name of gene with LOSS or AMP respectively
      this_gene = rownames(output_data[grepl(this_gene, rownames(output_data)),])[maybe_cnv]
      message(paste0("Found hotspot for gene ", this_gene, " that also has CNV data, processing ..."))
      output_data[,(output_data[c(this_gene),]>0 & output_data[c(hot),]==1)][c(this_gene, hot),][c(this_gene),] = 0
    }else{ # otherwise just use the gene nae
      message(paste0("Found hotspot for gene ", this_gene, ", processing ..."))
      output_data[,(output_data[c(this_gene),]>0 & output_data[c(hot),]==1)][c(this_gene, hot),][c(this_gene),] = 0
    }
    # if the above statement work, then there should be no overlaps between hotspot and any other mutations
    # for the same gene
    if(length(output_data[,(output_data[c(this_gene),]>0 & output_data[c(hot),]==1)][c(this_gene, hot),][c(this_gene),])==0){
      message("Success")
    }else{
      message(paste0("Problem occured with the ", feat_with_hotspot_data, " and the gene ", this_gene, "and there is still overlap between mutation and hotspot."))
      break
    }
  }

  # did user provide any features they would like to drop from matrix?
  if(!missing(drop_these_features)){
    message("You provided features to be dropped from matrix, removing them ...")
    output_data =
      output_data[!rownames(output_data) %in% drop_these_features,]
    message("Success")
  }

  # drop features that are occuring at a very low frequency
  low_feat = which(rowSums(output_data) <= floor(ncol(output_data)*min_feature_percent))
  if (length(low_feat)>0){
    message(paste0 ("There are ", length(low_feat), " features underrepresented and not meeting the minimum frequency of ", min_feature_percent))
    print(names(low_feat))
    output_data = output_data[-c(low_feat),]
  }else{
    message(paste0 ("There are ", length(low_feat), " features not meeting the minimum frequency of ", min_feature_percent))
    message("Proceeding without dropping any feature ...")
  }

  # are there any samples with 0 features? Yes, 1 exome and 1 genome
  samples_with_zero_feat = which(colSums(output_data) == 0)
  if (length(samples_with_zero_feat)>0){
    message(paste0 ("There are ", length(samples_with_zero_feat), " samples with no features and they will be dropped from matrix: "))
    print(names(samples_with_zero_feat))
    output_data = output_data[, -c(samples_with_zero_feat)]
  }else{
    message("All samples in your matrix are having at least one feature. Proceeding without dropping any samples ...")
  }





  # convert to matrix explicitly to make it NMF-input compatible
  output_data = as.matrix(output_data)

  return(output_data)

}


#' Standartize the chr prefix in a vector of chromosome names based on projection.
#' Internal helper function not meant for out-of-package use.
#'
#'
#' @param incoming_vector Input vector of any length with chromosome names.
#' @param projection Projection to which chr prefix should be standardized.
#'
#' @return A vector of chromosome names with prefix standardized to projection
#'
#' @examples
#' these_chrs = c(8, "13", "chr4", "chrY")
#' standardize_chr_prefix(these_chrs, projection = "hg38")
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


#' Calculate proportion of genome altered by CNV.
#'
#' `calculate_pga` returns a data.frame with estimated proportion of genome altered for each sample.
#'
#' This function calculates the percent of genome altered (PGA) by CNV. It takes into account the total length of
#' sample's CNV and relates it to the total genome length to return the proportion affected by CNV. The input is expected to be seg file.
#' The path to a local SEG file can be provided instead. If The custom seg file is provided, the minimum required columns are
#' sample, chrom, start, end, and log.ratio. The function can work with either individual or multi-sample seg file. The telomeres are always
#' excluded from calculation, and centromeres/sex chromosomes can be optionally included or excluded.
#'
#' @param this_seg Input data frame of seg file.
#' @param seg_path Optionally, specify the path to a local seg file.
#' @param projection Argument specifying the projection of seg file, which will determine chr prefix, chromosome coordinates, and genome size. Default is grch37, but hg38 is also accepted.
#' @param cutoff The minimum log.ratio for the segment to be considered as CNV. Default is 0.56, which is 1 copy. This value is expected to be positive float of log.ratio for both deletions and amplifications.
#' @param exclude_sex Boolean argument specifying whether to exclude sex chromosomes from calculation. Default is TRUE.
#' @param exclude_centromeres Boolean argument specifyng whether to exclude centromeres from calculation. Default is TRUE.
#'
#' @return A data frame of sample_id and a float in the range [0..1] indicating the fraction of genome altered by CNV.
#' @export
#' @import data.table tidyverse
#'
#' @examples
#' sample_seg = get_sample_cn_segments(this_sample_id = "14-36022T") %>% rename("sample"="ID")
#' calculate_pga(this_seg = sample_seg)
#' calculate_pga(this_seg = sample_seg, exclude_sex = FALSE)
#'
#' multi_sample_seg = rbind(get_sample_cn_segments(this_sample_id = "14-36022T"),
#'                          get_sample_cn_segments(this_sample_id = "BLGSP-71-21-00243-01A-11E")) %>%
#'                          rename("sample"="ID")
#' GAMBLR::calculate_pga(this_seg = multi_sample_seg)
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
    this_seg = read_tsv(seg_path)
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

#' Adjust ploidy for samples with CNV data.
#'
#' `adjust_ploidy` returns a seg file with log.ratios adjusted to the overall sample ploidy.
#'
#' This function adjusts ploidy of the sample using the percent of genome altered (PGA). The PGA is calculated internally, but can also be optionally provided as data frame
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
#' @export
#' @import tidyverse
#'
#' @examples
#' sample_seg = get_sample_cn_segments(this_sample_id = "14-36022T") %>% rename("sample"="ID")
#' adjust_ploidy(this_seg = sample_seg)
#'
#' multi_sample_seg = rbind(get_sample_cn_segments(this_sample_id = "14-36022T"),
#'                          get_sample_cn_segments(this_sample_id = "BLGSP-71-21-00243-01A-11E")) %>%
#'                          rename("sample"="ID")
#' adjust_ploidy(this_seg=multi_sample_seg)
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
    this_seg = read_tsv(seg_path)
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
    if (!"sample" %in% colnames(pga)) {
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


#' Helper function called by fancy_multisample_ideo, for sub-setting copy number information based on segments avaialble in cn data
#'
#' @param cn_segments DF with copy number segments, usually retrieved from get_sample_cn_segments.
#' @param include_2 Optional parameter for including or ommit CN state == 2.
#' @param samplen Numeric value that annotates the sample order.
#'
#' @return Nothing.
#'
#' @examples
#' cn_states = get_sample_cn_segments(multiple_samples = TRUE, sample_list = c("00-15201_tumorA", "HTMCP-01-06-00422-01A-01D"), streamlined = FALSE)
#' subset_cnstates(cn_segments = cn_states, samplen = 1)
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


#' Compare segmented data for multiple samples.
#'
#' `cnvKompare` returns a list in variable data formats allowing to evaluate concordance of CNV data between multiple samples.
#'
#' This function will compare CNV data between samples with multiple time points. It can also handle same-sample comparison
#' between different CNV callers if sample ID is specified in unique fashion. For groups with more than 2 samples,
#' optionally the pairwise comparisons can be performed. The comparison is made based on the internally calculated score,
#' which reflects percentage of each cytoband covered by CNV (rounded to the nearest 5%) and its absolute CN. Optionally,
#' the heatmap of cnvKompare scores can be returned. In addition, the function will return all concordant and discordant cytobands.
#' Finally, the time series plot of CNV log ratios will be returned for all lymphoma genes, with further functionality to subset
#' it to a panel of genes of interest.
#'
#' @param patient_id Specify patient_id to retrieve sample ids from GAMBL metadata.
#' @param sample_ids Optionally, specify sample ids for comparison.
#' @param this_seg Optional input data frame of seg file. Must adhere to seg format.
#' @param seg_path Optionally, specify the path to a local seg file. Must adhere to seg format.
#' @param genes_of_interest Provide specific genes to be displayed on the time-series plot.
#' @param projection Argument specifying the projection of seg file, which will determine coordinates of the cytobands. Default is grch37, but hg38 is also accepted.
#' @param ignore_cytoband_labels Cytobands to be ignored. By default, "acen", "gvar", "stalk" are excluded.
#' @param max_overlap For time-series plot, how many maximum overlapping points are allowed?
#' @param min_concordance Integer value from 0 to 100 to indicate the minimum required similarity between cytobands to be considered concordant. The dafult is 90 (90%).
#' @param pga_cutoff Minimum PGA for the sample to adjust ploidy. Default is 0.05 (5%).
#' @param exclude_sex Boolean argument specifying whether to exclude sex chromosomes from calculation. Default is FALSE.
#' @param return_heatmap Boolean argument specifying whether to return a heatmap of cnvKompare scores. Default is TRUE.
#' @param compare_pairwise Boolean argument specifying whether to perform pairwise comparisons is there are more than 2 time points in the group. Default is TRUE.
#'
#' @return A list of overall and pairwise percent concordance, concordant and discordant cytobands, comparison heatmap of cnvKompare socres, and time series ggplot object.
#' @export
#' @import tidyverse data.table circlize ComplexHeatmap ggrepel
#'
#' @examples
#' cnvKompare(patient_id = "00-14595", genes_of_interest = c("EZH2", "TP53", "MYC", "CREBBP", "GNA13"))
#' cnvKompare(patient_id = "13-26835", genes_of_interest = c("EZH2", "TP53", "MYC", "CREBBP", "GNA13"), projection = "hg38")
cnvKompare = function(patient_id,
                      sample_ids,
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
  if (missing(patient_id) & missing(sample_ids)) {
    stop("Please provide patient id or sample ids for comparison.")
  }

  # retrieve sample ids if only patient id is specified
  if (missing(sample_ids)) {
    sample_ids = get_gambl_metadata()
    sample_ids = dplyr::filter(sample_ids, patient_id == {{ patient_id }})
    sample_ids = pull(sample_ids, sample_id)
    message(paste0(
      "Found ",
      length(sample_ids),
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
    these_samples_seg = read_tsv(seg_path) %>%
      `names<-`(c(ID, chrom, start, end, LOH_flag, log.ratio)) %>%
      dplyr::mutate(CN = (2 * 2 ^ log.ratio))
  } else if (!missing(this_seg)) {
    these_samples_seg = this_seg %>%
      `names<-`(c(ID, chrom, start, end, LOH_flag, log.ratio)) %>%
      dplyr::mutate(CN = (2 * 2 ^ log.ratio))
  } else {
    message("Retreiving the CNV data using GAMBLR ...")
    these_samples_seg = get_sample_cn_segments(multiple_samples = TRUE,
                                               sample_list = sample_ids,
                                               from_flatfile = TRUE,
                                               projection = projection,
                                               with_chr_prefix = TRUE)
  }

  these_samples_seg = these_samples_seg  %>%
    dplyr::filter(ID %in% sample_ids) %>% # if user-provided seg, ensure only samples of comparison are present
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
  if (compare_pairwise & length(sample_ids) > 2) {
    message("Performing pairwise comparisons ...")

    # generate all possible combinations
    possible_combinations = apply(combn(sample_ids, 2), 2, paste, collapse =
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


#' Transform input maf columns to allow for usage of dplyr verbs
#'
#' @param maf_df input MAF data frame.
#'
#' @return maf_df with transofmred columns
#' @export
#' @import tidyverse
#'
#' @examples
#' ssm_sample = get_ssm_by_sample(this_sample_id = "HTMCP-01-06-00485-01A-01D", tool_name = "slims-3", projection = "grch37")
#' clean_maf = cleanup_maf(maf_df = ssm_sample)
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

#' Complement maf with missing samples.
#'
#' @param incoming_maf The initial MAF data frame to be supplemented with missing samples.
#' @param these_samples_metadata The metadata data frame that contains Tumor_Sample_Barcode column with ids to be present in the complemented MAF.
#'
#' @return maf_df with complemented Tumor_Sample_Barcode and other columns ready to be used downstream
#' @export
#'
#' @examples
#' small_maf = get_coding_ssm(limit_cohort = "dlbcl_reddy", seq_type = "capture") %>% dplyr::filter(Hugo_Symbol=="MYC")
#' reddy_meta = get_gambl_metadata(seq_type_filter = "capture") %>% dplyr::filter(cohort=="dlbcl_reddy")
#' complete_maf = supplement_maf(incoming_maf = small_maf, these_samples_metadata = reddy_meta)
#'
supplement_maf <- function(incoming_maf,
                           these_samples_metadata){
  missing_sample_ids = setdiff(these_samples_metadata$Tumor_Sample_Barcode,
                               incoming_maf$Tumor_Sample_Barcode)
  missing_sample_maf = incoming_maf %>%
    dplyr::filter(Tumor_Sample_Barcode == "Imaginary Sample ID") %>%
    add_row(Tumor_Sample_Barcode = missing_sample_ids,
           Hugo_Symbol = "GARBAGE",
           Chromosome = ifelse(stringr::str_detect(incoming_maf$Chromosome[1],
                                                   "chr"),
                               "chr1",
                               "1"),
           Start_Position = 1,
           End_Position = 1,
           Variant_Classification = "Missense_Mutation"
           )
  full_maf = rbind(incoming_maf,
                   missing_sample_maf)
  return(full_maf)
}
