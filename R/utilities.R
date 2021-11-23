
#' Get a MAF that is just the variants unique to one of two flavours of variant calls available
#'
#' @param these_sample_ids
#' @param flavour1
#' @param flavour2
#'
#' @return
#' @export
#'
#' @examples
compare_mutation_flavour = function(these_sample_ids,flavour1="clustered",flavour2=""){
  these_dfs = list()
  for(this_sample_id in these_sample_ids){
    message(this_sample_id)
    maf1 = get_ssm_by_sample(this_sample_id,flavour=flavour1)
    maf2  = get_ssm_by_sample(this_sample_id,flavour=flavour2)
    maf1_only = intersect_maf(maf1,maf2)
    these_dfs[[this_sample_id]]=maf1_only
  }

  this_maf = rbindlist(these_dfs,use.names = TRUE)
  return(this_maf)
}

#' perform set operations on two MAFs
#'
#' @param maf1
#' @param maf2
#' @param set_returned
#'
#' @return
#' @export
#'
#' @examples
intersect_maf = function(maf1,maf2,set_returned="maf1_only"){
  if(set_returned=="maf1_only"){
    maf_set = dplyr::filter(maf1,!Start_Position %in% maf2$Start_Position)
  }else if(set_returned == "maf2_only"){
    maf_set = dplyr::filter(maf2,!Start_Position %in% maf1$Start_Position)
  }
  return(maf_set)
}

#' Tabulate mutation status for non-silent SSMs for a set of genes
#'
#' @param gene_symbols List of gene symbols for which the mutation status will be tabulated. If not provided, lymphoma genes will be returned by default.
#' @param these_samples_metadata The matedata for samples of interest to be included in the returned matrix. Only the column "sample_id" is required. If not provided, the matrix is tabulated for all available samples as default.
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations.
#' @param include_hotspots Logical parameter indicating whether hotspots object should also be tabulated. Default is TRUE.
#' @param from_flatfile Integer value indicating minimal recurrence level
#' @param review_hotspots Logical parameter indicating whether hotspots object should be reviewed to include functionally relevant mutations or rare lymphoma-related genes. Default is TRUE.
#' @param ... Other parameters accepted by the review_hotspots() function
#'
#' @return
#' @export
#'
#' @examples
#' coding_tabulated_df = get_coding_ssm_status(gene_symbols=c("MYC","KMT2D"))
#' coding_tabulated_df = get_coding_ssm_status() #all lymphoma genes from bundled NHL gene list
get_coding_ssm_status = function(gene_symbols,
                                  these_samples_metadata,
                                  from_flatfile=TRUE,
                                  include_hotspots=TRUE,
                                  recurrence_min = 5,
                                  review_hotspots=TRUE,
                                  genes_of_interest = c("FOXO1", "MYD88", "CREBBP"),
                                  genome_build = "hg19"){
  if(missing(gene_symbols)){
    message("defaulting to all lymphoma genes")
    gene_symbols = pull(lymphoma_genes,Gene)
  }
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata()
  }

  # call it once so the object can be reused if user wants to annotate hotspots
  coding_ssm = get_coding_ssm(from_flatfile=from_flatfile)

  coding = coding_ssm %>%
    dplyr::filter(Hugo_Symbol %in% gene_symbols &
                    Variant_Classification != "Synonymous") %>%
    dplyr::select(Tumor_Sample_Barcode,Hugo_Symbol) %>%
    dplyr::rename("sample_id"="Tumor_Sample_Barcode","gene"="Hugo_Symbol") %>%
    unique() %>%
    mutate(mutated=1)
  samples_table = dplyr::select(these_samples_metadata,sample_id)
  wide_coding = pivot_wider(coding,names_from = "gene",
                            values_from="mutated",values_fill = 0)
  #complete(wide_coding,fill=list("sample_id"=samples_table$sample_id))
  all_tabulated = left_join(samples_table,wide_coding)
  all_tabulated = all_tabulated %>% replace(is.na(.), 0)

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
              dplyr::select(Tumor_Sample_Barcode,Hugo_Symbol, hot_spot) %>%
              dplyr::rename("sample_id"="Tumor_Sample_Barcode","gene"="Hugo_Symbol") %>%
              dplyr::mutate(gene=paste0(gene, "HOTSPOT")) %>%
              unique() %>%
              dplyr::mutate(mutated=ifelse(hot_spot=="TRUE", 1, 0)) %>%
              dplyr::filter(mutated==1) %>%
              dplyr::select(-hot_spot)

    # long to wide hotspots, samples are tabulated with 0 if no hotspot is detected
    wide_hotspots = pivot_wider(hotspots,names_from = "gene",
                          values_from="mutated",values_fill = 0)
    # join with the ssm object
    all_tabulated = left_join(all_tabulated,wide_hotspots)
    all_tabulated = all_tabulated %>% replace(is.na(.), 0)
    # make SSM and hotspots non-redundant by giving priority to hotspot feature and setting SSM to 0
    for (hotspot_site in colnames(wide_hotspots)[grepl("HOTSPOT", colnames(wide_hotspots))]){
      message(hotspot_site)
          this_gene = gsub("HOTSPOT", "", hotspot_site)
          redundant_features = all_tabulated %>% dplyr::select(starts_with(this_gene))
          # if not both the gene and the hotspot are present, go to the next iteration
          if(ncol(redundant_features)!=2) next
          message("OK")
          # if both gene and it's hotspot are in the matrix, give priority to hotspot feature
          all_tabulated[(all_tabulated[,this_gene]>0 & all_tabulated[,paste0(this_gene, "HOTSPOT")]==1),][,c(this_gene, paste0(this_gene, "HOTSPOT"))][,this_gene] = 0
    }

  }

  return(all_tabulated)
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
trim_scale_expression <- function(x){
  quants <- unname(quantile(x, probs = c(0.05, 0.95),na.rm=TRUE))
  x <- ifelse(x < quants[1], quants[1], x)
  x <- ifelse(x > quants[2], quants[2], x)
  x <- (x - quants[1]) / (quants[2] - quants[1])
  return(x)
}

#' Count hypermutated bins and generate heatmaps/cluster the data
#'
#' @param regions Vector of regions in the format "chr:start-end"
#' @param region_df Data frame of regions with four columns (chrom,start,end,gene_name)
#' @param slide_by How far to shift before starting the next window
#' @param window_size The width of your sliding window
#' @param min_count_per_bin
#' @param min_bin_recurrence How many samples a bin must be mutated in to retain in the visualization
#' @param min_bin_patient How many bins must a patient mutated in to retain in the visualization
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process (or full metadata)
#' @param region_padding How many bases will be added on the left and right of the regions to ensure any small regions are sufficiently covered by bins
#' @param metadataColumns What metadata will be shown in the visualization
#' @param sortByColumns Which of the metadata to sort on for the heatmap
#' @param cluster_rows_heatmap Optional parameter to enable/disable clustering of each dimension of the heatmap
#' @param cluster_cols_heatmap
#' @param customColour Optional named list of named vectors for specifying all colours for metadata. Can be generated with map_metadata_to_colours
#' @param show_gene_colours Optional logical argument indicating whether regions should have associated colours plotted as annotation track of heatmap
#' @param legend_row Fiddle with these to widen or narrow your legend
#' @param legend_col Fiddle with these to widen or narrow your legend
#' @param legend_col Accepts one of "horizontal" (default) or "vertical" to indicate in which direction the legend will be drawn
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details)
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#'
#' @return
#' @export
#'
#' @examples
get_mutation_frequency_bin_matrix = function(regions,
                                  regions_df,
                                  these_samples_metadata,
                                  region_padding= 1000,
                                  metadataColumns=c("pathology"),
                                  sortByColumns=c("pathology"),
                                  expressionColumns=c(),
                                  orientation = "sample_rows",
                                  customColour = NULL,
                                  slide_by=100,
                                  window_size=500,
                                  min_count_per_bin = 3,
                                  min_bin_recurrence = 5,
                                  min_bin_patient = 0,
                                  region_fontsize=8,
                                  cluster_rows_heatmap = FALSE,
                                  cluster_cols_heatmap = FALSE,
                                  show_gene_colours=FALSE,
                                  legend_row=3,
                                  legend_col=3,
                                  legend_direction="horizontal",
                                  from_indexed_flatfile=FALSE,
                                  mode="slms-3"){

    if(missing(regions)){
      if(missing(regions_df)){
        regions_df = grch37_ashm_regions #drop MYC and BCL2
        regions_df = grch37_ashm_regions %>%
          dplyr::filter(!gene %in% c("MYC","BCL2","IGLL5"))
      }
      regions = unlist(apply(regions_df,1,function(x){paste0(x[1],":",as.numeric(x[2])-region_padding,"-",as.numeric(x[3])+region_padding)})) #add some buffer around each
    }
    dfs = lapply(regions,function(x){calc_mutation_frequency_sliding_windows(
    this_region=x,drop_unmutated = TRUE,
    slide_by = slide_by,plot_type="none",window_size=window_size,
    min_count_per_bin=min_count_per_bin,return_count = TRUE,
    metadata = these_samples_metadata,
    from_indexed_flatfile=from_indexed_flatfile, mode=mode)})

  all= do.call("rbind",dfs)
  #add a fake bin for one gene and make every patient not mutated in it (to fill gaps)
  fake = these_samples_metadata %>% dplyr::select(sample_id) %>% mutate(bin="1_chrN") %>% mutate(mutated=0)
  all = bind_rows(all,fake)
  completed = complete(all,sample_id,bin,fill = list(mutated = 0))
  widened = pivot_wider(completed,names_from=sample_id,values_from=mutated)
  widened_df = column_to_rownames(widened,var="bin")

  #meta_show = metadata %>% select(sample_id,pathology,lymphgen) %>%
  if(length(expressionColumns)>0){
    these_samples_metadata = these_samples_metadata %>%
      mutate(across(all_of(expressionColumns), ~ trim_scale_expression(.x)))
  }
  meta_show = these_samples_metadata %>% select(sample_id,all_of(metadataColumns)) %>%
    arrange(across(all_of(sortByColumns))) %>%
    dplyr::filter(sample_id %in% colnames(widened_df)) %>%
    column_to_rownames(var="sample_id")
  message(paste("starting with",length(colnames(widened_df)),"patients"))
  patients_show = colnames(widened_df)[which(colSums(widened_df)>=min_bin_patient)]
  message(paste("returning matrix with",length(patients_show),"patients"))
  meta_show = dplyr::filter(meta_show,rownames(meta_show) %in% patients_show)
  to_show = widened_df[which(rowSums(widened_df)>min_bin_recurrence),patients_show]

  bin_col_fun = colorRamp2(c(0, 3, 6, 9),
                       c("white", "orange","red","purple"))
  to_show_t = t(to_show)
  meta_show_t = meta_show[rownames(to_show_t),]
  lg_cols = get_gambl_colours("lymphgen")
  path_fun = function(x){
    path_cols = get_gambl_colours("pathology")
    lg_cols = get_gambl_colours("lymphgen")

    return(unname(path_cols[x]))
  }
  path_cols = get_gambl_colours("pathology")

  # assign bins back to regions for better annotation

  assign_bins_to_region = function(bin_names,rdf){
    bin_df = data.frame(bin_name=bin_names)

    separated = bin_df %>%
      separate(bin_name,into=c("start","chrom")) %>%
      mutate(start = as.integer(start)) %>%
      mutate(end=start+1)

    separated$bin_name = bin_names
    colnames(rdf)[c(1:3)]=c("chrom","start","end")
    rdf = mutate(rdf,start=start-1500) %>% mutate(end=end+1500)
    regions.dt = as.data.table(rdf)

    setkey(regions.dt,chrom,start,end)
    bin.dt = as.data.table(separated)
    setkey(bin.dt,chrom,start,end)
    bin_overlapped = foverlaps(bin.dt,regions.dt,mult="first") %>%
      as.data.frame() %>% select(bin_name,gene) %>% column_to_rownames(var="bin_name")
    return(bin_overlapped)
  }
  #regions_df = grch37_ashm_regions
  if(is.null(customColour)){
    meta_cols = map_metadata_to_colours(metadataColumns,these_samples_metadata = meta_show,as_vector = F)

  }else{
    meta_cols = customColour
  }
  col_fun=circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  for(exp in expressionColumns){
    meta_cols[[exp]] = col_fun
  }
  bin_annot = assign_bins_to_region(bin_names=colnames(to_show_t),rdf=regions_df)
  heatmap_legend_param = list(title = "Bin value",nrow=legend_row, ncol=legend_row,
                         legend_direction = legend_direction)

  annotation_legend_param = list(nrow=legend_row,
                            ncol=legend_col,
                            direction=legend_direction)

  if(orientation == "sample_rows"){
    row_annot = HeatmapAnnotation(df=meta_show,show_legend = T,
                                  which = 'row',
                                  col=meta_cols,
                                  annotation_legend_param = annotation_legend_param)
    if(show_gene_colours){
        col_annot = HeatmapAnnotation(df=bin_annot,show_legend = F,
                                  which = 'col',
                                  annotation_legend_param = annotation_legend_param)
    }else{
        col_annot = HeatmapAnnotation(value=anno_empty(border = FALSE))
    }
    Heatmap(to_show_t[rownames(meta_show),rownames(bin_annot)],
           cluster_columns = cluster_cols_heatmap,
           cluster_rows=cluster_rows_heatmap,
           col=bin_col_fun,
           bottom_annotation = col_annot,
           left_annotation = row_annot,
           show_row_names = F,show_column_names = F,
           column_split = factor(bin_annot$gene),
           #row_split = factor(meta_show$pathology),
           column_title_gp = gpar(fontsize=region_fontsize),
           column_title_rot = 90,
           row_title_gp = gpar(fontsize=10),
           heatmap_legend_param = heatmap_legend_param)
  }else{
    col_annot = HeatmapAnnotation(df=meta_show,show_legend = T,
                                  which = 'col',
                                  col=meta_cols,
                                  annotation_legend_param = annotation_legend_param)
    if(show_gene_colours){
      row_annot = HeatmapAnnotation(df=bin_annot,show_legend = F,
                                  which = 'row',
                                  annotation_legend_param = annotation_legend_param)
    }else{
      row_annot = rowAnnotation(value=anno_empty(border = FALSE))
    }

    Heatmap(to_show[rownames(bin_annot),rownames(meta_show)],show_heatmap_legend = F,
            cluster_columns = cluster_rows_heatmap,
            cluster_rows=cluster_cols_heatmap,
            col=bin_col_fun,
            bottom_annotation = col_annot,
            left_annotation = row_annot,
            show_row_names = F,show_column_names = F,
            row_split = factor(bin_annot$gene),
            #row_split = factor(meta_show$pathology),
            row_title_gp = gpar(fontsize=region_fontsize),
            row_title_rot = 0,
            column_title_gp = gpar(fontsize=8),
            heatmap_legend_param = heatmap_legend_param)
  }


}

#' Count the number of mutations in a sliding window across a region for all samples. Unlikely to be used directly in most cases. See get_mutation_frequency_bin_matrix instead
#'
#' @param chromosome
#' @param start_pos
#' @param end_pos
#' @param metadata
#' @param return_format
#' @param classification_column Only used for plotting
#' @param plot_type Set to true for a plot of your bins. By default no plots are made.
#' @param min_count_per_bin
#' @param return_count
#' @param drop_unmutated This may not currently work properly.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details)
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#' @return
#' @export
#' @import dplyr data.table ggplot2 cowplot
#'
#' @examples

calc_mutation_frequency_sliding_windows =
  function(this_region,chromosome,start_pos,end_pos,
           metadata,slide_by=100,
           window_size=1000,
           plot_type = "none",
           return_format="long-simple",
           min_count_per_bin=3,
           return_count = FALSE,
           drop_unmutated=FALSE,
           classification_column="lymphgen",
           from_indexed_flatfile=FALSE,
           mode="slms-3"){


  max_region = 1000000
  if(missing(metadata)){
    metadata = collate_results(join_with_full_metadata = TRUE)
  }
  if(missing(this_region)){
    this_region =paste0(chromosome,":",start_pos,"-",end_pos)
  }else{
    chunks = region_to_chunks(this_region)
    #print(chunks)
    chromosome = chunks$chromosome
    start_pos=as.numeric(chunks$start)
    end_pos=as.numeric(chunks$end)
  }
  region_size = end_pos - start_pos
  if(region_size < max_region){
    message(paste("processing bins of size",window_size,"across",region_size,"bp region"))
  }else{
    message(paste("CAUTION!\n",region_size,"exceeds maximum size recommended by this function."))

  }
  windows = data.frame(start=seq(start_pos,end_pos,by=slide_by)) %>%
    mutate(end=start+window_size-1)


  #use foverlaps to assign mutations to bins
  windows.dt = as.data.table(windows)


  region_ssm = GAMBLR::get_ssm_by_region(region=this_region,streamlined = TRUE, from_indexed_flatfile=from_indexed_flatfile, mode=mode) %>%
    dplyr::rename(c("start"="Start_Position","sample_id"="Tumor_Sample_Barcode")) %>%
    mutate(mutated=1)

  region.dt = region_ssm %>%
    dplyr::mutate(start=as.numeric(as.character(start)),
                  end=start+1,
                  end=as.numeric(as.character(end))) %>%
    dplyr::relocate(start, .before=end) %>%
    as.data.table()
  setkey(windows.dt,start,end)
  setkey(region.dt,start,end)

  windows_overlap = foverlaps(windows.dt,region.dt) %>%
    dplyr::filter(!is.na(start)) %>%
    dplyr::rename(c("window_start"="i.start","mutation_position"="start")) %>%
    dplyr::select(-i.end,-end,-mutation_position) %>% as.data.frame()

  windows_tallied_full = windows_overlap %>%
    group_by(sample_id,window_start) %>%
    tally() %>%
    dplyr::filter(n>=min_count_per_bin) %>%
    arrange(sample_id) %>% as.data.frame()
  windows_tallied = windows_tallied_full


  all_samples = pull(metadata,sample_id) %>% unique()
  num_samples = length(all_samples)
  lg_cols = get_gambl_colours("lymphgen")
  path_cols = get_gambl_colours("pathology")
  annos = data.frame(window_start = rep(start_pos,num_samples),
                     sample_id=factor(all_samples))
  annos = left_join(annos,metadata,by="sample_id")
  windows_tallied = left_join(metadata,windows_tallied,by="sample_id")
  windows_tallied$classification = factor(windows_tallied[,classification_column],levels=unique(windows_tallied[,classification_column]))
  if(drop_unmutated){
    windows_tallied = windows_tallied %>% dplyr::filter(!is.na(n))
  }
  if(classification_column== "lymphgen"){
    windows_tallied = arrange(windows_tallied,pathology,lymphgen)
    annos = arrange(annos,pathology,lymphgen)
  }else{
    windows_tallied = arrange(windows_tallied,classification)
    annos = arrange(annos,classification)

  }
  annos$sample_id = factor(annos$sample_id,levels=unique(annos$sample_id))
  windows_tallied$sample_id = factor(windows_tallied$sample_id,levels=unique(windows_tallied$sample_id))

  if(plot_type %in% c("points","point")){
    #p = ggplot2::ggplot(windows_tallied,aes(x=window_start,y=sample_id,colour=bcl2_ba)) +
    #geom_point(alpha=0.5) + theme(axis.text=element_text(size=4))
    #add a bin at position 1 for pathology
    windows_tallied = dplyr::filter(windows_tallied,!is.na(window_start))

    p = ggplot2::ggplot() +
      geom_point(data=annos,aes(x=window_start,y=sample_id,
                                colour=pathology)) +
      geom_point(data=windows_tallied,aes(x=window_start,y=sample_id,
                  colour=classification)) +
      theme(axis.text=element_text(size=4)) +
      scale_colour_manual(values=c(lg_cols,path_cols))

  }else if(plot_type == "tile"){
    p = windows_tallied %>%
      ggplot(aes(x=window_start,y=sample_id,fill=n)) +
      geom_tile() +scale_fill_gradient(low = "orange", high = "red", na.value = NA) +
      theme_cowplot()
  }
  if(plot_type != "none"){
    print(p)
  }
  if(return_count){
    windows_tallied = mutate(windows_tallied,bin = paste0(window_start,"_",chromosome)) %>%
      mutate(mutated=n)
  }else{
    #binary mutated or not
    windows_tallied = mutate(windows_tallied,bin = paste0(window_start,"_",chromosome)) %>%
    mutate(mutated=1)
  }
  a = dplyr::select(windows_tallied,sample_id,bin,mutated)
  completed = complete(a,sample_id,bin,fill = list(mutated = 0))
  widened = pivot_wider(completed,names_from=sample_id,values_from=mutated)
  if(return_format == "long"){
    return(windows_tallied)
  }else if(return_format == "long-simple"){
    #just return the columns needed for completing and making a wide matrix for many regions
    windows_simple = dplyr::select(windows_tallied,sample_id,bin,mutated)
    return(windows_simple)

  }else{
    return(widened)
  }
}

#' Write bedpe format data frame to a file that will work with IGV
#'
#' @param sv_df data frame of bedpe formatted SV data
#' @param filename File name (will be written to results/icgc_dart/misc/FILENAME )
#' @param add_chr_prefix Whether to force chr to be added to chromosome names
#'
#' @return
#' @export
#'
#' @examples
sv_to_bedpe_file = function(sv_df,filename="something.bedpe",add_chr_prefix=TRUE){
  #add chr prefix if missing
  if(add_chr_prefix){
    if(!grepl("chr",region_sv$CHROM_A[1])){
      sv_df = mutate(sv_df,CHROM_A=paste0("chr",CHROM_A)) %>%
        mutate(CHROM_B=paste0("chr",CHROM_B))
    }
  }
  bed_file="/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/misc/"
  bed_file = paste0(bed_file,filename)
  write.table(sv_df,file=bed_file,sep="\t",quote=F,row.names=F,col.names=F)
}

#' Parse a region string into chomosome, start and end
#'
#' @param region A region string e.g. "chrX:12345-678910"
#'
#' @return A named list
#' @export
#'
#' @examples
#' chr_start_end = region_to_chunks("chr1:1111-2222")
region_to_chunks = function(region){
  region = unname(region)
  region = gsub(",","",region)
  #format is chr6:37060224-37151701
  split_chunks = unlist(strsplit(region,":"))
  chromosome = split_chunks[1]
  startend = unlist(strsplit(split_chunks[2],"-"))
  qstart=startend[1]
  qend=startend[2]
  return(list(chromosome=chromosome,start=qstart,end=qend))
}

#' Write an oncomatrix from a MAF File for further plotting. This is meant to be run by individuals who have access to data sets to
#' "sanitize" a subset of data for subsequent use by them or others who don't have permission to access the raw data.
#' Example: User J has full permissions for ICGC data and has read permissions on a MAF file. User B needs to make some oncoplots
#' and/or perform some statistical analysis on the frequency and assortment of mutations in that data set but doesn't need all the details.
#' User J can run this function on a maf file and provide the path of the output to user B.
#'
#' @param mutation_maf_path Provide either the full path to a MAF file or
#' @param mutation_maf_data Otherwise provide a data frame of the MAF data
#' @param output_oncomatrix Optionally provide the path for your sanitized output file (otherwise it writes to working directory)
#' @param genes_keep Specify which genes you want to remain in the output
#' @param genes_drop Optionally specify which genes to drop (this doesn't mean all other genes will remain. Maftools decides that part)
#'
#' @return The full path to the oncomatrix file (a matrix with Variant_Classification or Multi_Hit indicating coding mutation status per patient)
#' @export
#'
#' @examples
#' lymph_genes = lymphoma_genes$Gene #note, this will be used by default if the user wants to be lazy
#' secure_maf = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/slms-3_vcf2maf_current/level_3/final_merged_grch37.CDS.maf"
#' safe_oncomatrix_path = sanitize_maf_data(mutation_maf_path=secure_maf,genes_keep=lymph_genes)
#'
sanitize_maf_data = function(mutation_maf_path,mutation_maf_data,output_oncomatrix,genes_keep,genes_drop=c()){
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

  maftools::oncoplot(maf_o,genes=genes_keep,writeMatrix = T,removeNonMutated = F)  #writes to working directory
  if(!missing(output_oncomatrix)){
    #rename it
    file.dplyr::rename("onco_matrix.txt",output_oncomatrix)
  }else{
    output_oncomatrix=paste0(getwd(),"/onco_matrix.tsv")
  }
  message(paste("your data is in:",output_oncomatrix))
  return(output_oncomatrix)
}

#' Annotate MAF-like data frome with a hot_spot column indicating recurrent mutations
#'
#' @param mutation_maf A data frame in MAF format
#' @param analysis_base Base name for hot spot output directory
#'
#' @return The same data frame with one additional column "hot_spot"
#' @export
#'
#' @examples
#' hot_ssms = annotate_hotspots(all_ssm)
#' hot_maf = read.maf(hot_ssms)
#' oncoplot(hot_maf,genes=c("MEF2B","TP53","MYD88"),additionalFeature = c("hot_spot",TRUE))
annotate_hotspots = function(mutation_maf,recurrence_min = 5,analysis_base=c("FL--DLBCL","BL--DLBCL"),p_thresh=0.05){
  hotspot_info = list()
  for(abase in analysis_base){
    base_path="/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/icgc_dart/oncodriveclustl-0.0/99-outputs/webversion/"
    clust_full_path = paste0(base_path,abase,"/NO_SILENT_MUTS/",abase,"_clusters_results.tsv")
    all_full_path = paste0(base_path,abase,"/NO_SILENT_MUTS/",abase,"_elements_results.txt")
    clust_hotspot=read_tsv(clust_full_path)
    all_hotspot=read_tsv(all_full_path)
  clustered_hotspots = clust_hotspot %>% dplyr::select(-RANK) %>% dplyr::filter(N_SAMPLES>recurrence_min & P < p_thresh)
  #clustered_hotspots = clust_hotspot %>% dplyr::select(-RANK)
  #Use the max and min coordinate to annotate all non-silent mutations in that region
  #all_ssm %>% filter(Chromosome=="19" & Start_Position>= 19260033 & Start_Position <= 19260045)

    arranged = clustered_hotspots %>% separate_rows(COORDINATES,convert=TRUE) %>%
      group_by(SYMBOL,MAX_COORD) %>% arrange(COORDINATES)

    mins = arranged %>% slice_head() %>% dplyr::rename("START"="COORDINATES")
    maxs = arranged %>% slice_tail() %>% dplyr::rename("END"="COORDINATES")
    hotspot_ranges = left_join(mins, dplyr::select(maxs, c(MAX_COORD,END)), by = c("SYMBOL","MAX_COORD"))
    hotspot_info[[abase]]=hotspot_ranges
  }
  merged_hotspot = do.call("rbind",hotspot_info)  %>% ungroup()

  long_hotspot = merged_hotspot %>% dplyr::select(MAX_COORD,CHROMOSOME,START,END) %>%
   pivot_longer(c(START,END),names_to="which",values_to="COORDINATE") %>% dplyr::select(-which)
  #again take highest and lowest value for each MAX_COORD
  starts = long_hotspot %>% group_by(MAX_COORD) %>% arrange(COORDINATE) %>% slice_head()
  ends = long_hotspot %>% group_by(MAX_COORD) %>% arrange(COORDINATE) %>% slice_tail()
  long_hotspot = bind_rows(starts,ends)
  filled_coords = long_hotspot  %>% group_by(MAX_COORD) %>% arrange(MAX_COORD,COORDINATE) %>%
  complete(COORDINATE = seq(COORDINATE[1], COORDINATE[2]))  %>%
    fill(CHROMOSOME, .direction = "up") %>% dplyr::rename("Start_Position"="COORDINATE") %>%
    dplyr::rename("Chromosome"="CHROMOSOME") %>% ungroup()
  filled_coords = mutate(filled_coords,hot_spot=TRUE)
  #just the ssms that match these coordinates!
  hot_ssms = left_join(mutation_maf,filled_coords,by=c("Chromosome","Start_Position"))
  return(hot_ssms)
}

#' Annotate MAF-like data frome with a hot_spot column indicating recurrent mutations
#'
#' @param annotated_maf A data frame in MAF format that has hotspots annotated using function annotate_hotspots().
#' @param genes_of_interest List of genes for hotspot review. Currently only FOXO1, MYD88, and CREBBP are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF file. The default is hg19 genome build.
#'
#' @return The same data frame with reviewed column "hot_spot"
#' @export
#' @import dplyr
#'
#' @examples
#' hot_ssms = review_hotspots(annotate_hotspots(get_coding_ssm()), genes_of_interest=c("CREBBP"))

review_hotspots = function(annotated_maf, genes_of_interest=c("FOXO1", "MYD88", "CREBBP"), genome_build="hg19"){

  # check genome build because CREBBP coordinates are hg19-based or hg38-based
  coordinates <- list()
  if (genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
    coordinates$start <- 3785000
    coordinates$end <- 3791000
  }else if(genome_build %in% c("hg38", "grch38", "GRCh38")){
    coordinates$start <- 3734999
    coordinates$end <- 3740999
  }else{
    stop("The genome build specified is not currently supported. Please provide MAF file in one of the following cordinates: hg19, grch37, hs37d5, GRCh37, hg38, grch38, or GRCh38")
  }

  # check that at least one of the currently supported genes are present
  if (sum(c("FOXO1", "MYD88", "CREBBP") %in% genes_of_interest)<1){
      stop("Currently only FOXO1, MYD88, and CREBBP are supported. Please specify one of these genes.")
  }

  # notify user that there is limited number of genes currently supported
  if (sum(c("FOXO1", "MYD88", "CREBBP") %in% genes_of_interest)>1 & length(genes_of_interest) > 3 ){
      print("Currently only FOXO1, MYD88, and CREBBP are supported. By default only these genes from the supplied list will be reviewed.")
  }

  if("FOXO1" %in% genes_of_interest){
      annotated_maf <- annotated_maf %>%
          dplyr::mutate(hot_spot=ifelse(Hugo_Symbol=="FOXO1" & HGVSp_Short == "p.M1?", "TRUE" , hot_spot))
  }
  if("CREBBP" %in% genes_of_interest){
      annotated_maf <- annotated_maf %>%
          dplyr::mutate(hot_spot=ifelse(Hugo_Symbol=="CREBBP" & Start_Position > coordinates$start & End_Position < coordinates$end & Variant_Classification == "Missense_Mutation", "TRUE" , hot_spot))
  }
  if("MYD88" %in% genes_of_interest){
      annotated_maf <- annotated_maf %>%
          dplyr::mutate(hot_spot=ifelse(Hugo_Symbol=="MYD88" & HGVSp_Short %in% c("p.L273P", "p.L265P"), "TRUE" , hot_spot))
  }
  return(annotated_maf)
}


#' Make a UCSC-ready custom track file from SV data
#
#' @param sv_bedpe A bedpe formatted data frame of SVs
#' @param output_file A bed file with UCSC custom header
#'
#' @return nothing
#' @export
#'
#' @examples
#' all_sv = get_manta_sv()
#' sv_to_custom_track(all_sv,output_file="GAMBL_sv_custom_track.bed")
sv_to_custom_track = function(sv_bedpe,output_file,is_annotated=TRUE,sv_name="all"){
  #browser position chr7:127471196-127495720
  #browser hide all
  #track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"
  #chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
  if(is_annotated){
  #reduce to a bed-like format
  sv_data1 = mutate(sv_bedpe,annotation=paste0(chrom1,":",start1,"_",fusion)) %>%
    dplyr::select(chrom2,start2,end2,tumour_sample_id,annotation,fusion)
  sv_data2 = mutate(sv_bedpe,annotation=paste0(chrom2,":",start2,"_",fusion)) %>%
    dplyr::select(chrom1,start1,end1,tumour_sample_id,annotation,fusion)
  print(head(sv_data1))
  print(head(sv_data2))
  colnames(sv_data1)=c("chrom","start","end","sample_id","annotation","fusion")
  colnames(sv_data2)=c("chrom","start","end","sample_id","annotation","fusion")
  sv_data = bind_rows(sv_data1,sv_data2)
  sv_data = mutate(sv_data,end=end+10)
  }else{
    sv_data = mutate(sv_bedpe,annotation=paste0(1,":",2)) %>%
      dplyr::select(1,2,6,tumour_sample_id,annotation)
    #sv_data2 = mutate(sv_bedpe,annotation=paste0(4,":",5)) %>%
    #  dplyr::select(1,2,3,tumour_sample_id,annotation)
    colnames(sv_data)=c("chrom","start","end","sample_id","annotation")
    #colnames(sv_data2)=c("chrom","start","end","sample_id","annotation")
  }

  if(!grepl("chr",sv_data[,1])){
    #add chr
    sv_data[,1] = unlist(lapply(sv_data[,1],function(x){paste0("chr",x)}))
  }
  coo_cols=get_gambl_colours("COO")
  path_cols = get_gambl_colours("pathology")
  all_cols=c(coo_cols,path_cols)
  colour_df = data.frame(coo=names(all_cols),colour=all_cols)
  rgb_df = data.frame(t(col2rgb(all_cols))) %>%
    mutate(consensus_coo_dhitsig=names(all_cols)) %>%
    unite(col="rgb",red,green,blue,sep = ",")
  meta=get_gambl_metadata() %>% dplyr::select(sample_id,"consensus_coo_dhitsig",pathology) %>%
    mutate(consensus_coo_dhitsig=if_else(consensus_coo_dhitsig=="NA",pathology,consensus_coo_dhitsig))

  samples_coloured = left_join(meta, rgb_df)
  sv_bed_coloured = left_join(sv_data,samples_coloured) %>% arrange(pathology)

  #chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
  write_bed = function(coloured_svs,sv_name,output_file_base){
    data_bed = coloured_svs %>% mutate(details=paste0(annotation,"_",sample_id)) %>%
      mutate(score=0,strand="+",end=end+1,start1=start,end1=end) %>%
      dplyr::select(chrom, start,end, details,score,strand,start1,end1,rgb) %>%
      dplyr::filter(!is.na(rgb)) %>% unique()
    header_content=paste0('track name="GAMBL SVs ',sv_name, '" description="SV breakpoints ', sv_name, '" visibility=2 itemRgb="On"\n')
    cat(header_content,file=output_file)
    message(paste("writing to",output_file))
    tabular = write.table(data_bed,file=output_file,quote=F,sep="\t",row.names=F,col.names = F,append=TRUE)
  }
  write_bed(sv_bed_coloured,sv_name=sv_name)

}

#' Convert a maf-formatted data frame into a bed custom track file for UCSC
#'
#' @param maf_data Either a maf loaded from disk or from the database using a get_ssm function
#' @param output_file Name for your new bed file that can be uploaded as a custom track to UCSC
#'
#' @return Nothing
#' @export
#'
#' @examples
#' maf_to_custom_track(my_maf_data,"/home/rmorin/private/some_mutations.bed")
maf_to_custom_track = function(maf_data,output_file){
  #browser position chr7:127471196-127495720
  #browser hide all
  #track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"
  #chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0

  #reduce to a bed-like format
  maf_data = dplyr::select(maf_data,Chromosome,Start_Position,End_Position,Tumor_Sample_Barcode)
  colnames(maf_data)=c("chrom","start","end","sample_id")
  if(!grepl("chr",maf_data[,1])){
    #add chr
    maf_data[,1] = unlist(lapply(maf_data[,1],function(x){paste0("chr",x)}))
  }
  lymphgen_cols=get_gambl_colours("lymphgen")
  colour_df = data.frame(lymphgen=names(lymphgen_cols),colour=lymphgen_cols)
  rgb_df = data.frame(t(col2rgb(lymphgen_cols))) %>%
    mutate(lymphgen=names(lymphgen_cols)) %>%
    unite(col="rgb",red,green,blue,sep = ",")
  meta=get_gambl_metadata() %>% dplyr::select(sample_id,lymphgen)
  samples_coloured = left_join(meta, rgb_df)
  maf_bed = maf_data %>% mutate(score=0,strand="+",end=end+1,start1=start,end1=end)
  maf_coloured = left_join(maf_bed,samples_coloured,by="sample_id") %>% dplyr::select(-lymphgen) %>% dplyr::filter(!is.na(rgb))
  cat('track name="GAMBL mutations" description="Mutations from GAMBL" visibility=2 itemRgb="On"\n',file=output_file)
  tabular = write.table(maf_coloured,file=output_file,quote=F,sep="\t",row.names=F,col.names = F,append=TRUE)
}


#' Bring together all derived sample-level results from many GAMBL pipelines
#'
#' @return A table keyed on biopsy_id that contains a bunch of per-sample results from GAMBL
#' @export
#' @import tidyverse config
#'
#' @examples
#' everything_collated = collate_results(join_with_full_metadata=TRUE)
collate_results = function(sample_table,
                           write_to_file=FALSE,
                           join_with_full_metadata = FALSE,
                           case_set,sbs_manipulation="",seq_type_filter="genome"){
  # important: if you are collating results from anything but WGS (e.g RNA-seq libraries) be sure to use biopsy ID as the key in your join
  # the sample_id should probably not even be in this file if we want this to be biopsy-centric
  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter = seq_type_filter) %>%
      dplyr::select(sample_id,patient_id,biopsy_id)
  }
  message("this will be slow until collate_ssm_results is modified to cache its result")
  #edit this function and add a new function to load any additional results into the main summary table
  sample_table = collate_ssm_results(sample_table=sample_table)
  sample_table = collate_sv_results(sample_table=sample_table)
  sample_table = collate_curated_sv_results(sample_table=sample_table)
  sample_table = collate_ashm_results(sample_table=sample_table)
  sample_table = collate_nfkbiz_results(sample_table=sample_table)
  sample_table = collate_csr_results(sample_table=sample_table)
  sample_table = collate_ancestry(sample_table=sample_table)
  sample_table = collate_sbs_results(sample_table=sample_table,sbs_manipulation=sbs_manipulation)
  sample_table = collate_derived_results(sample_table=sample_table)

  if(write_to_file){
    output_file = config::get("table_flatfiles")$derived
    output_base = config::get("repo_base")
    output_file = paste0(output_base,output_file)
    write_tsv(sample_table,file=output_file)
  }
  #convenience columns bringing together related information
  if(join_with_full_metadata){
  full_meta = get_gambl_metadata(seq_type_filter = seq_type_filter)
  full_table = left_join(full_meta,sample_table)
  full_table = full_table %>% mutate("MYC_SV_any"=case_when(ashm_MYC > 3 ~ "POS",
                                            manta_MYC_sv == "POS" ~ "POS",
                                            ICGC_MYC_sv == "POS" ~ "POS",
                                            myc_ba == "POS" ~ "POS",
                                            TRUE ~ "NEG"))
  full_table = full_table %>% mutate("BCL2_SV_any"=case_when(ashm_BCL2 > 3 ~ "POS",
                                                            manta_BCL2_sv == "POS" ~ "POS",
                                                            ICGC_BCL2_sv == "POS" ~ "POS",
                                                            bcl2_ba == "POS" ~ "POS",
                                                            TRUE ~ "NEG"))
  full_table =full_table %>% mutate("DoubleHitBCL2"=ifelse(BCL2_SV_any == "POS" & MYC_SV_any=="POS","Yes","No"))
  return(full_table)
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
#' gambl_results_derived = collate_derived_results(samples_df)
collate_derived_results = function(sample_table,from_flatfile=FALSE){

  if(from_flatfile){
    message("not implemented YET")
  }else{
    database_name = config::get("database_name")

    con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    derived_tbl = dplyr::tbl(con,"derived_data") %>% as.data.frame()
  }
  derived_tbl = derived_tbl %>%
    dplyr::select(where( ~!all(is.na(.x)))) %>% as.data.frame() #drop the columns that are completely empty
  print(derived_tbl)
  sample_table = dplyr::left_join(sample_table,derived_tbl)
  print(sample_table)
  return(sample_table)
}


#' Collate a few CSR annotations, including MiXCR
#'
#' @param sample_table A data frame with sample_id as the first column
#'
#' @return The sample table with additional columns
#' @export
#' @import tidyverse
#'
#' @examples
#' gambl_results_derived = collate_csr_results(gambl_results_derived)
collate_csr_results = function(sample_table){
   csr = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-nthomas/results/icgc_dart/mixcr_current/level_3/mixcr_genome_CSR_results.tsv"))
   sm_join = inner_join(sample_table,csr,by=c("sample_id"="sample"))
   pt_join = inner_join(sample_table,csr,by=c("patient_id"="sample"))
   complete_join <- bind_rows(pt_join, sm_join) %>%
     bind_rows(dplyr::filter(sample_table, !patient_id %in% c(pt_join$patient_id, sm_join$patient_id))) %>% unique()
  return(complete_join)
}

#' Compute some summary statistics based on SSM calls
#'
#' @param sample_table
#'
#' @return
#' @export
#'
#' @examples
collate_ssm_results = function(sample_table,from_flatfile=TRUE){
  coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del",
                   "In_Frame_Ins","Missense_Mutation","Nonsense_Mutation",
                   "Nonstop_Mutation","Splice_Region",
                   "Splice_Site","Targeted_Region","Translation_Start_Site")
  #iterate over every sample and compute some summary stats from its MAF
  if(from_flatfile){
    base_path = config::get("project_base")
    #test if we have permissions for the full gambl + icgc merge
    maf_partial_path = config::get("results_filatfiles")$ssm$all$full
    maf_path = paste0(base_path,maf_partial_path)
    maf_permissions = file.access(maf_path,4)
    if(maf_permissions == -1){
      #currently this will only return non-ICGC results
      maf_partial_path = config::get("results_filatfiles")$ssm$gambl$full
      base_path = config::get("project_base")
      #default is non-ICGC
      maf_path = paste0(base_path,maf_partial_path)
    }
    muts=fread_maf(maf_path)
    mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
    message(paste("mutations from",mutated_samples,"samples"))
  }
  #get tally of total per sample
  muts = muts %>% dplyr::rename("sample_id"="Tumor_Sample_Barcode")
  muts = mutate(muts,vaf=t_alt_count/(t_alt_count+t_ref_count))
  muts_count = dplyr::select(muts,sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("total_ssm"="n")
  sample_table = left_join(sample_table,muts_count)
  muts_mean = muts %>% dplyr::select(sample_id,vaf) %>%
    group_by(sample_id) %>%
    summarize(mean_vaf=mean(vaf))
  coding_mut = dplyr::filter(muts,Variant_Classification %in% coding_class)
  coding_mut_count = coding_mut %>%
    dplyr::select(sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("coding_ssm"="n")

  sample_table = left_join(sample_table,muts_mean)
  sample_table = left_join(sample_table,coding_mut_count)
  #check for coding SSMs in lymphoma genes
  coding_nhl = coding_mut %>%
    dplyr::filter(Hugo_Symbol %in% lymphoma_genes$Gene)
  coding_nhl_count = coding_nhl %>% group_by(sample_id) %>% tally() %>%
    dplyr::rename("driver_ssm"="n")
  return(sample_table)
}



#' Collate all SV calls from the genome data and summarize for main oncogenes of interest per sample
#'
#' @param sample_table A data frame with sample_id as the first column
#'
#' @return The sample table with additional columns
#' @export
#' @import tidyverse
#'
#' @examples
#' gambl_results_derived = collate_curated_sv_results(gambl_results_derived)
collate_curated_sv_results = function(sample_table){
  path_to_files = config::get("derived_and_curated")
  project_base = config::get("project_base")
  #  "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
  manual_files = dir(paste0(project_base,path_to_files),pattern=".tsv")
  for(f in manual_files){
    full = paste0(project_base,path_to_files,f)
    this_data = suppressMessages(read_tsv(full,comment = "#"))
    #TO DO: fix this so it will join on biopsy_id or sample_id depending on which one is present
    sample_table = left_join(sample_table,this_data)
  }


  return(sample_table)
}

#' Annotate mutations with their copy number information
#'
#' @param this_sample Sample ID of the sample you want to annotate
#' @param seg_file_source Specify what copy number calling program the input seg file is from, as it handles ichorCNA differently than WisecondorX, Battenberg, etc.
#' @param coding_only Optional. set to TRUE to rescrict to only coding variants
#' @param from_flatfile Optional. Instead of the database, load the data from a local MAF and seg file
#' @param assume_diploid Optional. If no local seg file is provided, instead of defaulting to a GAMBL sample, this parameter annotates every mutation as copy neutral.
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
assign_cn_to_ssm = function(this_sample,
                            coding_only=FALSE,
                            from_flatfile=FALSE,
                            use_augmented_maf=FALSE,
                            maf_file,
                            seg_file,
                            seg_file_source="ichorCNA", 
                            assume_diploid=FALSE){

  database_name = config::get("database_name")
  project_base = config::get("project_base")
  #tool_name=config::get("analyses")$matched$copy_number


  #project_base = "/projects/nhl_meta_analysis_scratch/gambl/results_local/"
  coding_class = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")
  if(!missing(maf_file)){
    maf_sample = fread_maf(maf_file) %>%
      dplyr::mutate(Chromosome = gsub("chr","",Chromosome))
  }else if(from_flatfile){
    #get the genome_build for this sample
    bam_info = get_bams(this_sample)
    #message(paste("bams:",bam_info))
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
    message(paste("loading MAF:",this_sample_maf))
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
  #if(tool_name == "battenberg"){
  
  if(!missing(seg_file)){
    seg_sample = read_tsv(seg_file) %>%
      dplyr::mutate(size=end - start) %>%
      dplyr::filter(size > 100)
    colnames(seg_sample)[c(1:4)] = c("ID","chrom","start","end")
    seg_sample = seg_sample %>%
      dplyr::mutate(chrom = gsub("chr","",chrom)) %>%
      dplyr::rename(Chromosome=chrom,Start_Position=start,End_Position=end) %>%
      data.table::as.data.table()
    #print(seg_sample)
    data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    a = data.table::as.data.table(maf_sample)
  }else if(assume_diploid == TRUE){
    if(missing(seg_file)){
      print("WARNING: A seg file was not provided! Annotating all mutation calls as copy neutral")
    }
    a = data.table::as.data.table(maf_sample)
    a_diploid = dplyr::mutate(a, CN=2)
    return(a_diploid)
  }else if(from_flatfile){
      message(paste("fetching:",tool_name))
      battenberg_files = fetch_output_files(build=genome_build,base_path = "gambl/battenberg_current",tool="battenberg",search_pattern = ".igv.seg")

      battenberg_file = dplyr::filter(battenberg_files,tumour_sample_id==this_sample) %>%
        dplyr::pull(full_path) %>% as.character()
      message(paste("using flatfile:",battenberg_file))
      if(length(battenberg_file)>1){
        print("WARNING: more than one SEG found for this sample. This shouldn't happen!")
        battenberg_file = battenberg_file[1]
      }
      seg_sample = read_tsv(battenberg_file) %>%
        as.data.table() %>% dplyr::mutate(size=end - start) %>%
        dplyr::filter(size > 100) %>%
        dplyr::mutate(chrom = gsub("chr","",chrom)) %>%
        dplyr::rename(Chromosome=chrom,Start_Position=start,End_Position=end)
      
      data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
      a = data.table::as.data.table(maf_sample)
    }else{
      seg_sample = get_sample_cn_segments(sample_id=this_sample) %>%
        dplyr::mutate(size=end - start) %>%
        dplyr::filter(size > 100) %>%
        dplyr::mutate(chrom = gsub("chr","",chrom)) %>%
        dplyr::rename(Chromosome=chrom,Start_Position=start,End_Position=end) %>%
        data.table::as.data.table()
      
      data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
      a = data.table::as.data.table(maf_sample)
    
      #seg_table = config::get("results_tables")$copy_number
      #message(paste("querying",seg_table))
      #seg_sample = dplyr::tbl(con,seg_table) %>%
      #dplyr::filter(ID == this_sample) %>%
      #data.table::as.data.table() %>% dplyr::mutate(size=end - start) %>%
      #dplyr::filter(size > 100) %>%
      #dplyr::mutate(chrom = gsub("chr","",chrom)) %>%
      #dplyr::rename(Chromosome=chrom,Start_Position=start,End_Position=end)
    }
    #data.table::setkey(seg_sample, Chromosome, Start_Position, End_Position)
    
    if(!missing(seg_file_source)){
      if(seg_file_source=="ichorCNA"){
        #message("defaulting to ichorCNA format")
        seg_sample = dplyr::rename(seg_sample, c("log.ratio"="median", "CN"="copy.number"))
        a.seg = data.table::foverlaps(a, seg_sample, type="any")
        a$log.ratio = a.seg$log.ratio
        a$LOH = factor(a.seg$LOH_flag)
        a$CN = a.seg$CN
      }else{
        a.seg = data.table::foverlaps(a, seg_sample, type="any")
        a$log.ratio = a.seg$log.ratio
        a$LOH = factor(a.seg$LOH_flag)
        a = dplyr::mutate(a,CN=round(2*2^log.ratio))
        seg_sample = dplyr::mutate(seg_sample,CN=round(2*2^log.ratio))
        seg_sample$LOH_flag = factor(seg_sample$LOH_flag)
      }
    }
  
    
  #    mutate(a,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
  #    geom_point(aes(x=Start_Position,y=vaf,colour=CN),size=0.1)  +
  #    geom_segment(data=seg_sample,aes(x=Start_Position,xend=End_Position,y=CN,yend=CN,colour=LOH_flag)) +
  #    facet_wrap(~Chromosome,scales="free_x")
    if(!missing(seg_sample)){
      return(list(maf=a,seg=seg_sample))
    }

}


#' Estimate purity 
#'
#' @param in_maf Path to a local maf file
#' @param in_seg Path to a local corresponding seg file for the same sample ID as the input maf
#' @param seg_file_source Specify what copy number calling program the input seg file is from, as it handles ichorCNA differently than WisecondorX, battenberg, etc.
#' @param show_plots Optional. Show two faceted plots that display the VAF and purity distributions for each copy number state in the sample
#' @param assume_diploid Optional. If no local seg file is provided, instead of defaulting to a GAMBL sample, this parameter annotates every mutation as copy neutral.
#'
#' @return A list containing a data frame (MAF-like format) with the segmented absolute copy number data and three extra columns:
#' VAF is the variant allele frequency calculated from the t_ref_count and t_alt_count
#' Ploidy is the number of copies of an allele in the tumour cell
#' Final_purity is the finalized purity estimation per mutation after considering different copy number states and LOH events
#' @export
#' @import tidyverse
#'
#' @examples
#' tumour_purity = estimate_purity(in_maf="path/to/file.maf", in_seg="path/to/file.seg", seg_file_source="ichorCNA", show_plots=TRUE)
#' tumour_purity = estimate_purity(in_maf="path/to/file.maf", assume_diploid = TRUE)

estimate_purity = function(in_maf, 
                           in_seg, 
                           seg_file_source="ichorCNA",
                           show_plots=FALSE,
                           assume_diploid=FALSE
                           ){
  
  # Merge the CN info to the corresponding MAF file, uses GAMBLR function
  if(!missing(in_seg)){
    if(seg_file_source=="ichorCNA"){
      CN_new <- assign_cn_to_ssm(maf_file = in_maf, seg_file = in_seg, seg_file_source = "ichorCNA")$maf
    }else{
      CN_new <- assign_cn_to_ssm(maf_file = in_maf, seg_file = in_seg, seg_file_source = "")$maf
    }
  }else{
    # If no seg file was provided and assume_diploid paramtere is set to true, 
    if(assume_diploid){
      CN_new <- assign_cn_to_ssm(maf_file = in_maf, assume_diploid = TRUE)
    }
  }
  
  # Change any homozygous deletions (CN = 0) to 1 for calculation purposes
  CN_new$CN[CN_new$CN<1] = 1
  
  # Select only the relevant columns, add new columns for VAF, Ploidy, and Purity
  CN_new <- CN_new %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, t_ref_count, t_alt_count, CN) %>%
    dplyr::mutate(VAF = t_alt_count/(t_ref_count+t_alt_count)) %>%
    dplyr::mutate(Ploidy = ifelse(CN %in% 1:2, 1, 0)) %>% 
    dplyr::mutate(Purity = "") %>%
    tidyr::drop_na(CN) %>%
    tidyr::unite("Chrom_pos", Chromosome:Start_Position, sep=":", remove=FALSE)
  
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
    CN_max = CN_new %>% dplyr::filter(CN == i)
    if(i > 2){
      CN_max_dup = CN_max[rep(seq(nrow(CN_max)), CN_max$CN),]
      CN_max_dup$Ploidy = rep(seq(i),nrow(CN_max))
      CN_max_dup = dplyr::mutate(CN_max_dup, Purity = (CN*VAF)/Ploidy)
      indiv_CN[[as.character(i)]] = CN_max_dup
    }else(indiv_CN[[as.character(i)]] = CN_max)
  }
  
  # Merge all copy state tables together into one table
  merged_CN <- do.call("rbind", indiv_CN)
  
  # Filter for copy nnumber states 1 or 2
  # Calculate purity, and if the number is larger 1 (100%) use the VAF calue instead
  merged_CN_neut <- merged_CN %>%
    dplyr::filter(merged_CN$CN < 3) %>%
    dplyr::mutate(Purity = (CN*VAF)/Ploidy) %>%
    dplyr::mutate(Purity = ifelse(Purity > 1, VAF, Purity))
    
    # Calculate a temporary purity based on the mean of these purity values
    mean_neut_purity = mean(merged_CN_neut$Purity)
  
  # For CN of 3 or larger: 
    ## Filter for copy number states 1 or 2
    ## Calculate purity, and if the number is larger 1 (100%) use the VAF calue instead
    ## Group by chromosonal position, 
    ## For each group, choose only the purity value that is closest to the mean_neut_purity (the temporary purity calculated from the copy neutral values)
  if (is.na(mean_neut_purity)){
    merged_CN_gain <- merged_CN %>% 
      dplyr::filter(merged_CN$CN > 2) %>%
      dplyr::mutate(Purity = (VAF*2)) %>%
      dplyr::mutate(Purity = ifelse(Purity > 1, VAF, Purity)) %>%
      group_by(Chrom_pos) %>%
      slice_head()  %>%
      dplyr::mutate(Assumed_CN = 2) %>%
      dplyr::mutate(Assumed_ploidy = 1)
  }else{ 
    merged_CN_gain <- merged_CN %>% 
      dplyr::filter(merged_CN$CN > 2) %>%
      dplyr::mutate(Purity = (CN*VAF)/Ploidy) %>%
      dplyr::mutate(Purity = ifelse(Purity > 1, VAF, Purity)) %>%
      group_by(Chrom_pos) %>%
      slice(which.min(abs(Purity - mean_neut_purity)))  
  }

  # Bind both datatables back together (the first contains CNs 1 and/or 2, the second contains CNs 3 and/or higher)
  CN_final <- bind_rows(merged_CN_neut, merged_CN_gain)
  
  # Estimate the mean of all purity values from all available copy number states
  sample_purity_estimation = mean(CN_final$Purity)
  
  # Calculate CCF (cancer cell fraction) 2 ways:
   ## With the maximum purity estimation for the mutations in the same (largest value will be 1 for the mutation with the highest purity estimation)
   ## With the purity estimation of the sample in total. This will give values over 1 though which is problematic, as the maximum value should be 1 (100%) for clonal events
  CN_final <- CN_final %>%
    dplyr::mutate(CCF_mutation = Purity/max(Purity)) %>%
    dplyr::mutate(CCF_sample = Purity/sample_purity_estimation)
  
  output = list()
  
  if(show_plots){
    # Figure 1 : VAF distribution
    # Creates facet wraps showing the VAF distribution of CN-annotated mutations for each available copy number state
    VAF_plot <- CN_final %>% 
      ggplot(aes(x=VAF)) + 
      geom_histogram() + 
      facet_wrap(~CN) 
   
    # Figure 2: Final purity distribution
    Purity_plot <- CN_final %>% 
      ggplot(aes(x=Purity)) + 
      geom_histogram() + 
      facet_wrap(~CN)
     
     output[["VAF_plot"]] = VAF_plot
     output[["Purity_plot"]] = Purity_plot
  }
  
  output[["sample_purity_estimation"]] = sample_purity_estimation
  output[["CN_final"]] = CN_final
  
  return(output)
}


#' Refresh the contents of a database table
#'
#' @param table_name Name of table to refresh
#' @param connection Database connection object
#' @param file_path Path to the table contents to populate
#'
#' @return
#' @export
#'
#' @examples
#' refresh_full_table(table_x,con,file_x)
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
#' referesh_metadata_tables()
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
  tumour_samples = sample_df %>% dplyr::select(patient_id,sample_id,biopsy_id,seq_type,protocol) %>%
    dplyr::filter(!is.na(biopsy_id))
  n_samp_bio = tumour_samples %>% count() %>% pull(n)
  #2876 unique samples
  #check for any multiplicity of sample_id
  n_samp = tumour_samples %>% dplyr::select(-biopsy_id) %>% unique() %>% count() %>% pull(n)
  #should be the same number as above
  if(!n_samp == n_samp_bio){
    print("ERROR! some biopsies appear to have the same sample_id/protocol combination")
  }

  return(all_metadata_info)
}

collate_ancestry = function(sample_table,somalier_output){
  if(missing(somalier_output)){
    somalier_output="/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/gambl/somalier_current/02-ancestry/2020_08_07.somalier-ancestry.tsv"
  }
  somalier_all = read_tsv(somalier_output)
  somalier_all = mutate(somalier_all,sample_id = str_remove(`#sample_id`,pattern=":.+")) %>%
    dplyr::select(-`#sample_id`,-given_ancestry)
  somalier_all = dplyr::select(somalier_all,sample_id,predicted_ancestry,PC1,PC2,PC3,PC4,PC5)
  sample_table = left_join(sample_table,somalier_all)
  return(sample_table)
}


collate_extra_metadata= function(sample_table,file_path){
  file_path = "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/experiments/2021-04-21-Trios-MiXCR/trios_relapse_timing.tsv"
  extra_df = read_tsv(file_path)
  sample_table = left_join(sample_table,extra_df,by=c("sample_id"="biopsy_id"))
}

#' Bring in the results from mutational signature analysis
#'
#' @param sample_table A data frame with sample_id as the first column
#' @param file_path Optional path to SBS file
#'
#' @return A data frame with new columns added
#' @export
#' @import tidyverse
#'
#' @examples
#' collated = collate_sbs_results(sample_table=sample_table,sbs_manipulation=sbs_manipulation)
collate_sbs_results = function(sample_table,file_path,scale_vals=FALSE,sbs_manipulation=""){
  if(missing(file_path)){
    file_path = "/projects/rmorin_scratch/prasath_scratch/gambl/sigprofiler/gambl_hg38/02-extract/slms3.gambl.icgc.hg38.matched.unmatched/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt"
  }
  signatures = read.csv(file_path,sep="\t",header=1,row.names=1)
  rs=rowSums(signatures)
  cn=colnames(signatures)
  new_sig = signatures
  if(sbs_manipulation=="scale"){
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
    sbs = apply(signatures,2,function(x){x/rs}) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
  }
  else if(sbs_manipulation=="log"){
    #sbs1 = log(signatures[,"SBS1"]+1)
    #sbs5 = log(signatures[,"SBS5"]+1)
    #sbs9 = log(signatures[,"SBS9"]+1)
    #sbs8 = log(signatures[,"SBS8"]+1)
    #sbs12 = log(signatures[,"SBS12"]+1)
    #sbs17b = log(signatures[,"SBS17b"]+1)
    #sbs18 = log(signatures[,"SBS18"]+1)
    #sbs84 = log(signatures[,"SBS84"]+1)
    #sbs85 = log(signatures[,"SBS85"]+1)
    sbs = apply(signatures,2,function(x){log(x+1)}) %>%
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
    sbs = signatures %>% rownames_to_column("sample_id")
  }


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
#' sample_table = collate_nfkbiz_results(sample_table=sample_table)
collate_nfkbiz_results = function(sample_table){
  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>% dplyr::select(sample_id,patient_id,biopsy_id)
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
    sample_table = get_gambl_metadata() %>% dplyr::select(sample_id,patient_id,biopsy_id)
  }
  #just annotate BCL2, MYC and CCND1 hypermutation
  regions_df = data.frame(name=c("CCND1","BCL2","MYC"),
        region=c("chr11:69455000-69459900","chr18:60983000-60989000","chr8:128747615-128751834"))
  region_mafs = lapply(regions_df$region,function(x){get_ssm_by_region(region=x,streamlined = TRUE)})
  tibbled_data = tibble(region_mafs, region_name = regions_df$name)
  unnested_df = tibbled_data %>% unnest_longer(region_mafs)
  unlisted_df = mutate(unnested_df,start=region_mafs$Start_Position,sample_id=region_mafs$Tumor_Sample_Barcode) %>%
    dplyr::select(start,sample_id,region_name)
  tallied = unlisted_df %>% group_by(sample_id,region_name) %>%
    tally() %>%
    pivot_wider(values_from=n,names_from=region_name,values_fill=0,names_prefix="ashm_")
  #sample_table = mutate(sample_table,across(everything(), ~replace_na(.x, 0)))
  sample_table = left_join(sample_table,tallied,by="sample_id")


}

#' Determine and summarize which cases have specific oncogene SVs
#'
#' @param sample_table A data frame with sample_id as the first column
#' @param tool Name of tool (optional, default is manta)
#' @param oncogenes Which oncogenes to collate SVs from
#'
#' @return Data frame with additional columns ({tool}_{oncogene} and {tool}_{oncogene}_{partner})
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
    some_fusions = some_fusions %>% dplyr::select(tumour_sample_id,partner)  %>% mutate("{tool}_{oncogene_name}_partner" := partner) %>%
      dplyr::select(-partner)
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
#' @param classification (optionally request only colours for pathology, lymphgen, mutation or copy_number)
#'
#' @return A named vector of colour codes for lymphgen classes and pathology
#' @export
#' @import tidyverse
#'
#' @examples
#' lymphgen_cols=get_gambl_colours("lymphgen")
#' # be sure to install ggsci from https://github.com/morinlab/ggsci
#' # install_github("morinlab/ggsci")
get_gambl_colours = function(classification="all",alpha=1){
  all_colours = list()
  everything = c()
  blood_cols=ggsci::get_ash("blood")
  all_colours[["hmrn"]] <-c(
    "BCL2-MYC" = "#52000F",
    "BCL2"="#721F0F",
    "SOCS1/SGK1"="#D66B1F",
    "TET2/SGK1"="#C41230",
    "MYD88"="#3B5FAC",
    "NOTCH2"="#7F3293",
    "NOTCH1" = "#55B55E",
    "Other"="#ACADAF"
  )
  all_colours[["EBV"]] =  c("EBV-positive"="#7F055F",
                            "EBV-negative"="#E5A4CB",
                            "POS"="#7F055F",
                            "NEG"="#E5A4CB")
  all_colours[["BL"]] = c("M53-BL"="#A6CEE3",
                                         "DLBCL-1"="#721F0F",
                                         "IC-BL"="#45425A",
                                         "DGG-BL"="#33A02C",
                                         "DLBCL-2"="#FB9A99",
                                         "DLBCL-3"="#C41230")
  all_colours[["FL"]]=c(dFL="#99C1B9",cFL="#E8E46E")
  all_colours[["lymphgen"]] = c(
    "EZB-MYC" = "#52000F",
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
    "Other" = "#ACADAF"

  )
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

  all_colours[["pos_neg"]]=c(
    "POS"="#c41230",
    "NEG"="#E88873",
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
    "F"="#EF476F",
    "Female"="#EF476F")
  all_colours[["clinical"]]=ggsci::get_ash("clinical")
  all_colours[["pathology"]] = c(
      "B-ALL"="#C1C64B",
      "CLL"="#889BE5",
      "MCL"="#F37A20",
      "BL"="#926CAD",
      "mBL"="#34C7F4",
      "DLBCL-BL-like"="#34C7F4",
      "PMBL"= "#227C9D",
      "FL"="#EA8368",
      "COMFL"="#8BBC98",
      "DLBCL"="#479450",
      "HGBL-NOS"="#294936",
      "HGBL"="#294936",
      "HGBL-DH/TH"="#7A1616",
      "PBL" = "#E058C0",
      "CNS" = "#E2EF60",
      "THRLBCL" = "#A5F2B3",
      "MM"="#CC9A42",
      "SCBC"="#8c9c90",
      "UNSPECIFIED"="#cfba7c"
  )
  all_colours[["coo"]] = c(
    "ABC" = "#05ACEF",
    "UNCLASS" = "#05631E",
    "U" = "#05631E",
    "UNC" = "#05631E",
    "GCB"= "#F58F20",
    "DHITsig-"= "#F58F20",
    "DHITsigNeg"= "#F58F20",
    "DHITsig-IND" = "#003049",
    "DHITsig+" = "#D62828",
    "DHITsigPos" = "#D62828"
  )
  all_colours[["cohort"]] = c("Chapuy"="#8B0000",
                  "Arthur"= "#8845A8",
                  "Schmitz"= "#2C72B2")
  #print(all_colours)
  for(colslot in names(all_colours)){
    raw_cols=all_colours[[colslot]]
    raw_cols_rgb <- col2rgb(raw_cols)
    alpha_cols <- rgb(
      raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
      alpha = alpha * 255L, names = names(raw_cols),
      maxColorValue = 255L
    )
    names(alpha_cols)=names(raw_cols)
    all_colours[[colslot]]=alpha_cols
  }
  for(this_group in names(all_colours)){
    everything=c(everything,all_colours[[this_group]])
  }
  #return matching value from lowercase version of the argument if it exists
  lc_class = stringr::str_to_lower(classification)
  if(classification %in% names(all_colours)){
    return(all_colours[[classification]])
  }else if(lc_class %in% names(all_colours)){
    return(all_colours[[lc_class]])
  }else{
    #all=c(all_colours[["lymphgen"]],all_colours[["pathology"]],all_colours[["coo"]])
    #return(all)
    return(everything)
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


#' Using GISTIC2.0 outputs, perform Fisher's exact test to compare CNV frequencies between 2 groups
#'
#' @param gistic_lesions Path to the GISTIC2.0 all_lesions output file.
#' @param metadata Data frame containing sample ids and column with annotated data for the 2 groups of interest. All other columns are ignored. Currently, function exits if asked to compare more than 2 groups.
#' @param comparison Specify column annotating groups of interest.
#' @param fdr.method FDR method to adjust p values. Uses stats::p.adjust function, and therefore accepts its method for FDR ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). By default, this function uses "fdr".
#' @param fdr.cutoff Specify FDR significance cut-off. By default, this function uses 0.1.
#' @param text_size Size of the text on the forest plot of differentially enriched CNV.
#' @param blacklisted_regions Optionally, specify any descriptors (value from column `Descriptor` of GISTIC2.0 all_lesions output file) to filter out before amy comparisons are done. It is possible to specify list of multiple descriptors, for example, c("3p12.3", "12p13.2").
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
FtestCNV <- function(gistic_lesions, metadata, comparison, fdr.method="fdr", fdr.cutoff=0.1, text_size = 7, blacklisted_regions=NULL){
  # get groups of interest for comparison
  GROUPS.TO.COMPARE <- unique(metadata[,comparison])

  # quick check that only 2 groups are provided forr comparison
  if(length(GROUPS.TO.COMPARE) > 2){
    message("The current implementation of function only accepts 2 groups for comparison. You provided 3 groups.")
    message("Please modify metadata accordingly to compare only 2 groups.")
    message("Groups you provided are: ", paste(c(GROUPS.TO.COMPARE), collapse = ","))
    return(NULL)
  }

  # read lesions from gistic utput to collect event/case
  lesions <- readr::read_tsv(gistic_lesions, col_names = TRUE) %>%
    dplyr::filter(!grepl("CN", `Unique Name`)) %>%
    dplyr::select(-tail(names(.), 1), -`Residual q values after removing segments shared with higher peaks`, -`Broad or Focal`, -`Amplitude Threshold`, -`q values`) %>%
    dplyr::filter (! Descriptor %in% blacklisted_regions)

  # prepare metadata
  # save names of first few columns of gistic lesions file for future reference
  columns <- colnames(lesions)[1:5]
  # subset lesions file to only lesions metadata and samples of interest
  columns <- c(columns, pull(metadata, Tumor_Sample_Barcode))
  lesions <- subset(lesions, select=intersect(columns, colnames(lesions)))

  # subset ids for each category for future comparisons
  GROUP1 <- dplyr::pull(metadata %>% dplyr::filter(base::get(comparison)==GROUPS.TO.COMPARE[1]), Tumor_Sample_Barcode)
  GROUP2 <- dplyr::pull(metadata %>% dplyr::filter(base::get(comparison)==GROUPS.TO.COMPARE[2]), Tumor_Sample_Barcode)

  # subset lesions file for each group and count number of CNV
  GROUP1.CNV <- lesions[,colnames(lesions) %in% GROUP1] %>% dplyr::mutate(count=rowSums(.!=0))
  GROUP2.CNV <- lesions[,colnames(lesions) %in% GROUP2] %>% dplyr::mutate(count=rowSums(.!=0))

  ########## PART 1
  # compare significance of differences between 2 groups

  # first, get the counts of mutated samples for each group
  GROUP1_vs_GROUP2 <- as.data.frame(cbind(lesions[,1:5], GROUP1.CNV$count, GROUP2.CNV$count))
  colnames(GROUP1_vs_GROUP2)[6:7] <- c("Mutated_GROUP1", "Mutated_GROUP2")

  # second, get the counts of unmutated samples for each group
  GROUP1_vs_GROUP2 <- GROUP1_vs_GROUP2 %>% dplyr::mutate(Not_Mutated_GROUP1=length(GROUP1)-Mutated_GROUP1, .after = Mutated_GROUP1) %>%
    dplyr::mutate(Not_Mutated_GROUP2=length(GROUP2)-Mutated_GROUP2, .after = Mutated_GROUP2) %>%
    dplyr::mutate(Region=ifelse(grepl("Amplification", `Unique Name`), paste0(Descriptor, "_Amp"), paste0(Descriptor, "_Del")), .after = Descriptor)

  # row-by-row, calculate 2-tailed Fisher's exact test for each CNV
  GROUP1_vs_GROUP2 <- GROUP1_vs_GROUP2 %>%
    dplyr::filter(Mutated_GROUP1 / (Mutated_GROUP1 + Not_Mutated_GROUP1) > 0.05 | Mutated_GROUP2 / (Mutated_GROUP2 + Not_Mutated_GROUP2) > 0.05) %>%
    dplyr::mutate(pVal = apply(., 1, function(x) {tbl <- matrix(as.numeric(x[7:10]), ncol=2, byrow=T); fisher.test(tbl, alternative="two.sided")$p.value})) %>%
    dplyr::mutate(OddsRatio = apply(., 1, function(x) {tbl <- matrix(as.numeric(x[7:10]), ncol=2, byrow=T); log(fisher.test(tbl, alternative="two.sided")$estimate)})) %>%
    dplyr::mutate(LowConfInt = apply(., 1, function(x) {tbl <- matrix(as.numeric(x[7:10]), ncol=2, byrow=T); log(fisher.test(tbl, alternative="two.sided")$conf.int[1])})) %>%
    dplyr::mutate(HighConfInt = apply(., 1, function(x) {tbl <- matrix(as.numeric(x[7:10]), ncol=2, byrow=T); log(fisher.test(tbl, alternative="two.sided")$conf.int[2])}))

  # adjust FDR for multiple comparisons
  GROUP1_vs_GROUP2$FDR <- stats::p.adjust(GROUP1_vs_GROUP2$pVal, method = fdr.method)
  # filter for only those CNV that pass FDR cutoff
  GROUP1_vs_GROUP2.PASSED <- GROUP1_vs_GROUP2[GROUP1_vs_GROUP2$FDR<fdr.cutoff,]
  GROUP1_vs_GROUP2.PASSED <- GROUP1_vs_GROUP2.PASSED %>% dplyr::distinct(Region, .keep_all = TRUE)

  # rename columns to match with comparison groups and save resulting df as part of the output
  DISTINCT <- GROUP1_vs_GROUP2.PASSED %>% `colnames<-`(gsub("GROUP1", GROUPS.TO.COMPARE[1], colnames(.))) %>%
    `colnames<-`(gsub("GROUP2", GROUPS.TO.COMPARE[2], colnames(.)))

  message("Successfully completed step 1/3...")

  ########## PART 2
  # prepare data with CNV for more convenient downstream analyses, like survival

  # First, collect events for the group 1
  # this just extracts metadata for each CNV from lesions file
  REGIONS <- as.data.frame(lesions[,1:5]) %>%
    dplyr::mutate(Region=ifelse(grepl("Amplification", `Unique Name`), paste0(Descriptor, "_Amp"), paste0(Descriptor, "_Del")), .after = Descriptor)

  # Collect events for samples in group 1. Since same band can contain more than 1 peak, keep peak limits for future unique differentiation
  survival_object = 0
  # here, i is a CNV feature, anf j is a sample_id
  for(i in 1:length(REGIONS$Region)) {
    for(j in 1:1:ncol(GROUP1.CNV[,-ncol(GROUP1.CNV)])) {
      row = c(colnames(GROUP1.CNV[j]), REGIONS$Region[i], REGIONS$`Wide Peak Limits`[i], GROUP1.CNV[i,j])
      survival_object <- rbind(survival_object, row)
    }
  }

  # tidy output for convenience
  survival_object <- as.data.frame(survival_object)
  colnames(survival_object) <- c("sample_id", "Region", "Wide Peak Limits", "Event")
  survival_object <- survival_object[2:length(survival_object$sample),]
  survival_object <- survival_object %>% dplyr::filter(!Event==0) %>%
    mutate(Event = ifelse(grepl("Amp", Region), "GAIN", "LOSS"))

  # Create final object that will be concatenated with events for second group further down
  ALL.EVENTS <- survival_object

  # Collect events for samples in group 2
  survival_object = 0
  # here, i is a CNV feature, anf j is a sample_id
  for(i in 1:length(REGIONS$Region)) {
    for(j in 1:ncol(GROUP2.CNV[,-ncol(GROUP2.CNV)])) {
      row = c(colnames(GROUP2.CNV[j]), REGIONS$Region[i], REGIONS$`Wide Peak Limits`[i], GROUP2.CNV[i,j])
      survival_object <- rbind(survival_object, row)
    }
  }

  # tidy output for convenience
  survival_object <- as.data.frame(survival_object)
  colnames(survival_object) <- c("sample_id", "Region", "Wide Peak Limits", "Event")
  survival_object <- survival_object[2:length(survival_object$sample),]
  survival_object <- survival_object %>% dplyr::filter(!Event==0) %>%
    mutate(Event = ifelse(grepl("Amp", Region), "GAIN", "LOSS"))

  # Now, merge with the object prepared earlier and save it as part of the output
  CNV.EVENTS <- tidyr::unnest(rbind(ALL.EVENTS, survival_object))

  message("Successfully completed step 2/3...")

  ########## PART 3
  # visualize forest plot in a convenient way

  # calculate SE
  mergedPassed <- GROUP1_vs_GROUP2.PASSED %>%
    dplyr::mutate(SE = (HighConfInt - LowConfInt) / 5.95)
  # order in decreasing order for better visualization
  mergedPassed <- mergedPassed[order(mergedPassed$OddsRatio, decreasing = TRUE),]
  # prepare table that will be plotted
  study_table <- data.frame(
    name = mergedPassed[, "Region"],
    Events_GROUP1 = paste(mergedPassed$Mutated_GROUP1, "/", mergedPassed$Mutated_GROUP1+mergedPassed$Not_Mutated_GROUP1, sep = ""),
    Events_GROUP2 = paste(mergedPassed$Mutated_GROUP2, "/", mergedPassed$Mutated_GROUP2+mergedPassed$Not_Mutated_GROUP2, sep = ""))
  # rename columns to match with comparison groups
  study_table <- study_table %>% `colnames<-`(gsub("GROUP1", GROUPS.TO.COMPARE[1], colnames(.))) %>%
    `colnames<-`(gsub("GROUP2", GROUPS.TO.COMPARE[2], colnames(.)))

  # actual plot
  GRAPH <- metaviz::viz_forest(x = mergedPassed[, c("OddsRatio", "SE")], variant = "thick",
                      col = "Greys", xlab = "Log(OddsRatio)", annotate_CI = T,
                      type = "study_only",
                      study_table = study_table,
                      text_size = text_size,
                      table_headers = c("Region"))

  message("Successfully completed step 3/3...")

  OUTPUT <- list("DISTINCT"=DISTINCT, "CNV.EVENTS"=CNV.EVENTS, "GRAPH"=GRAPH)
  return(OUTPUT)
  message("Done!")
}



#' If some samples are missing from the matrix, add them with filled in 0 as value and normalize their ordering for consistency
#'
#' @param incoming_matrix A matrix or data frame that should be filled. Required parameter.
#' @param list_of_samples Vector specifying all desired samples to be present in the resulting matrix. Required parameter.
#' @param fill_in_values Value that will be used to fill in the matrix.
#' @param normalize_order Logical parameter specifying whether sample order should be according to the supplied list. Default is TRUE.
#' @param samples_in_rows Logical argument indicating whether samples are in rows or columns. Default assumes samples are in rows and columns are features.
#'
#' @return a data frame with maintained orientation (rows and columns) where samples from the supplied list are present and reordered according to the specified order
#' @export
#'
#' @examples
#' partial_matrix = get_coding_ssm_status(these_samples_metadata = (get_gambl_metadata(case_set = "BL--DLBCL") %>% filter(pairing_status=="unmatched")), include_hotspots = FALSE)
#' complete_matrix = complete_missing_from_matrix(partial_matrix, get_gambl_metadata() %>% pull(sample_id))
complete_missing_from_matrix = function(incoming_matrix,
                                        list_of_samples,
                                        fill_in_values = 0,
                                        normalize_order=TRUE,
                                        samples_in_rows=TRUE){

  # check for required arguments
  if (missing(incoming_matrix)){
      stop("Please provide initial matrix to fill.")
  }

  if (missing(list_of_samples)){
      stop("Please provide list of samples to complete the matrix and normalize order.")
  }

  # is samples are in columns, transpose the matrix so code below is generalizable
  if(!samples_in_rows){
    incoming_matrix = as.data.frame(incoming_matrix) %>% t()
  }

  matrix_with_all_samples <- rbind(incoming_matrix,
        matrix(fill_in_values:fill_in_values, # populate matrix with all 0
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
    matrix_with_all_samples = as.data.frame(matrix_with_all_samples) %>% t()
  }

  return(matrix_with_all_samples)
}
