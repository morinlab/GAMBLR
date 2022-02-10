
#' Generate a plot of all CN segments
#'
#' @param region
#' @param gene
#' @param these_samples_metadata
#' @param type
#' @param crop_segments
#' @param crop_distance
#'
#' @return
#' @export
#'
#' @examples
focal_cn_plot = function(region,gene,
                         these_samples_metadata,
                         type="gain",segment_size=1,
                         crop_segments = TRUE,
                         sort_by_annotation=c('pathology'),
                         crop_distance=100000000){
  if(!missing(gene)){
    region=gene_to_region(gene)
    chunks = region_to_chunks(region)
  }else{
    chunks = region_to_chunks(region)
  }
  if(type == "gain"){
    all_not_dip = get_cn_segments(region=region) %>%
      mutate(size=end-start) %>% dplyr::filter(CN>2)
  }else{
    all_not_dip = get_cn_segments(region=region) %>%
      mutate(size=end-start) %>% dplyr::filter(CN<2)
  }

  #crop start and end if they're further than crop_distance from your region
  all_not_dip = mutate(all_not_dip,left_distance=as.numeric(chunks$start)- start)
  all_not_dip = mutate(all_not_dip,right_distance=end - as.numeric(chunks$end))
  if(crop_segments){
    all_not_dip = mutate(all_not_dip,end=ifelse(right_distance >crop_distance,as.numeric(chunks$end)+crop_distance,end ))
    all_not_dip = mutate(all_not_dip,start=ifelse(left_distance >crop_distance,as.numeric(chunks$start)-crop_distance,start ))
  }
  #all_not_dip$ID=factor(all_not_dip$ID,levels=unique(all_not_dip$ID))
  all_not_dip = left_join(all_not_dip,these_samples_metadata,by=c("ID"="sample_id")) %>%
    dplyr::filter(!is.na(pathology))
  all_not_dip = all_not_dip %>% arrange(across(all_of(c(sort_by_annotation,"size"))))
  all_not_dip$ID=factor(all_not_dip$ID,levels=unique(all_not_dip$ID))
  ggplot(all_not_dip,aes(x=start,xend=end,y=ID,yend=ID,colour=lymphgen)) +
    geom_vline(aes(xintercept=as.numeric(chunks$start)),alpha=0.5,colour="red") +
    geom_segment(size=segment_size) + theme_cowplot() +
    theme(axis.text.y=element_blank())
}

#' Generate a more visually appealing and flexible lollipop plot
#'
#' @param maf_df A data frame containing the mutation data (from a MAF)
#' @param gene The gene symbol to plot
#' @param plot_title Optional (defaults to gene name)
#' @param plot_theme Options: default, blue, simple, cbioportal, nature, nature2, ggplot2, and dark
#'
#' @return
#' @export
#' @import g3viz
#'
#' @examples
pretty_lollipop_plot = function(maf_df,gene,plot_title,plot_theme="cbioportal"){
  if(missing(plot_title)){
    plot_title=gene
  }
  maf_df = maf_df %>% dplyr::filter(Hugo_Symbol==gene)
  #use the readMAF function (modified by Ryan) to parse/convert
  maf_df = g3viz::readMAF(maf.df=maf_df)
  chart.options <- g3Lollipop.theme(theme.name = plot_theme,title.text = plot_title)
  g3Lollipop(maf_df,
             gene.symbol = gene,
             plot.options = chart.options,
             output.filename = "default_theme")
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
                                             legendFontSize=10,
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
                              legend_direction = legend_direction,
                              labels_gp = gpar(fontsize = legendFontSize))

  annotation_legend_param = list(nrow=legend_row,
                                 ncol=legend_col,
                                 direction=legend_direction,
                                 labels_gp = gpar(fontsize = legendFontSize))

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


#' Plot a heatmap comparing the VAF of mutations in T1/T2 pairs
#'
#' @param maf1
#' @param maf2
#' @param vafcolname
#' @param patient_colname
#' @param these_samples_metadata
#' @param sortByColumns
#' @param metadata_columns
#' @param gene_orientation
#' @param annotate_zero
#' @param genes
#' @param top_n_genes
#' @param drop_unless_lowvaf
#' @param vaf_cutoff_to_drop
#' @param cluster_columns
#' @param cluster_rows
#'
#' @return
#' @export
#'
#' @examples
plot_mutation_dynamics_heatmap = function(maf1,maf2,vafcolname,
                                          patient_colname="patient_id",
                                          these_samples_metadata,
                                          sortByColumns,
                                          metadata_columns=c("sample_id"),
                                          gene_orientation="bottom",
                                          annotate_zero=FALSE,
                                          genes,
                                          top_n_genes,
                                          drop_unless_lowvaf=FALSE,
                                          vaf_cutoff_to_drop=0.05,
                                          cluster_columns=F,
                                          cluster_rows=F){
  if(missing(vafcolname)){
    t1_pair_mafdat = mutate(maf1,vaf=t_alt_count/(t_alt_count+t_ref_count))
    t2_pair_mafdat = mutate(maf2,vaf=t_alt_count/(t_alt_count+t_ref_count))
  }else{
    t1_pair_mafdat[,"vaf"]=t1_pair_mafdat[,vafcolname]
    t2_pair_mafdat[,"vaf"]=t2_pair_mafdat[,vafcolname]
  }
  if(!missing(sortByColumns)){
    these_samples_metadata = arrange(these_samples_metadata,across(sortByColumns))
  }
  #print(patient_colname)

  t1_pair_mafdat = t1_pair_mafdat %>% dplyr::rename("patient_id"=patient_colname)
  t2_pair_mafdat = t2_pair_mafdat %>% dplyr::rename("patient_id"=patient_colname)
  these_samples_metadata = these_samples_metadata %>% dplyr::rename("patient_id"=patient_colname)
  both_vaf_all = full_join(t1_pair_mafdat,t2_pair_mafdat,by=c("patient_id","Start_Position")) %>%
    dplyr::select(patient_id,HGVSp_Short.x,Hugo_Symbol.x,Hugo_Symbol.y,Tumor_Sample_Barcode.x,
                  Tumor_Sample_Barcode.y,HGVSp_Short.y,vaf.x,vaf.y,hot_spot.x,hot_spot.y) %>%
    mutate(ANNO=ifelse(is.na(vaf.x),HGVSp_Short.y,HGVSp_Short.x)) %>%
    mutate(GENE=ifelse(is.na(vaf.x),Hugo_Symbol.y,Hugo_Symbol.x)) %>%
    mutate(VAF1=ifelse(is.na(vaf.x),0,vaf.x)) %>%
    mutate(VAF2=ifelse(is.na(vaf.y),0,vaf.y)) %>%
    mutate(hot_spot.x=ifelse(is.na(hot_spot.x),0,1)) %>%
    mutate(hot_spot.y=ifelse(is.na(hot_spot.y),0,1)) %>%
    mutate(HOTSPOT=ifelse(hot_spot.x==1 | hot_spot.y==1,1,0)) %>%
    dplyr::filter(ANNO !="") %>%
    mutate(Mutation=paste(GENE,ANNO,sep="_")) %>%
    dplyr::select(patient_id,GENE,Mutation,VAF1,VAF2,ANNO,HOTSPOT)
  if(drop_unless_lowvaf){
    both_vaf_all = dplyr::filter(both_vaf_all,VAF1<vaf_cutoff_to_drop | VAF2<vaf_cutoff_to_drop)
  }
  if(!missing(genes)){
    both_vaf_all = dplyr::filter(both_vaf_all,GENE %in% genes)
  }
  both_vaf_all = mutate(both_vaf_all,unique_id=paste(patient_id,Mutation,sep="_")) %>%
    mutate(fold_change=log(VAF2+0.1)-log(VAF1+0.1)) %>%
    group_by(patient_id,GENE) %>%
    mutate(Number=paste(GENE,row_number(),sep="_"))
  print(head(both_vaf_all))
  h = both_vaf_all %>%
    select(patient_id,Number,fold_change) %>%
    pivot_wider(id_cols = patient_id,names_from=Number,values_from=fold_change) %>%
    column_to_rownames("patient_id")
  both_vaf_all = mutate(both_vaf_all,minVAF=ifelse(VAF1<VAF2,VAF1,VAF2))
  zeroes = both_vaf_all %>% select(patient_id,Number,minVAF) %>%
    pivot_wider(id_cols = patient_id,names_from=Number,values_from=minVAF) %>%
    column_to_rownames("patient_id")

  hotspots = both_vaf_all %>% select(patient_id,Number,HOTSPOT) %>%
    pivot_wider(id_cols = patient_id,names_from=Number,values_from=HOTSPOT) %>%
    column_to_rownames("patient_id")

  #zeroes = as.data.frame(zeroes)
  zeroes[is.na(zeroes)]=0.001
  hotspots[is.na(hotspots)]=-1
  hotspots[hotspots==0]=-1
  hotspots[hotspots==1]=0
  zeroes[zeroes>0]=1
  print(head(hotspots))

  #print(head(both_vaf_all))
  #print(head(these_samples_metadata))
  these_samples_metadata_rn =
    dplyr::filter(these_samples_metadata,patient_id %in% rownames(h)) %>%
    select(all_of(c("patient_id",metadata_columns))) %>%
    column_to_rownames("patient_id")

  la = HeatmapAnnotation(df = as.data.frame( these_samples_metadata_rn),
                         which="row")
  ta = HeatmapAnnotation(df = as.data.frame( these_samples_metadata_rn),
                         which="column")

  h[is.na(h)]=0
  cs=colSums(zeroes)
  ordered = names(cs[order(cs)])
  if(!missing(top_n_genes)){
    print(head(ordered,top_n_genes))
    genes=ordered[c(1:top_n_genes)]
  }
  col_fun = colorRamp2(c(0, 1), c("white", "red"))

  if(gene_orientation=="bottom"){
    if(annotate_zero){

    }else{
      H=Heatmap(h[rownames(these_samples_metadata_rn),],cluster_rows=F,
                cluster_columns=F,left_annotation = la)
    }
  }else{
    if(!missing(top_n_genes)){

      these_zeroes=t(zeroes[rownames(these_samples_metadata_rn),genes])
      these_zeroes=t(hotspots[rownames(these_samples_metadata_rn),genes])
      to_show = t(h[rownames(these_samples_metadata_rn),genes])

    }else{
      these_zeroes=t(zeroes[rownames(these_samples_metadata_rn),])
      these_zeroes=t(hotspots[rownames(these_samples_metadata_rn),])

      to_show = t(h[rownames(these_samples_metadata_rn),])
    }
    if(annotate_zero){


      H=Heatmap(to_show,

                layer_fun =
                  function(j, i, x, y, width, height, fill) {

                    v = pindex(these_zeroes, i, j)
                    l = v == 0
                    grid.points(x[l], y[l], pch = 16, size = unit(1, "mm"),
                                gp = gpar(col = "white"))
                    #grid.text(sprintf("%.1f", v[l]), x[l], y[l], gp = gpar(fontsize = 10))
                  },
                cluster_columns=cluster_columns,
                cluster_rows=cluster_rows,
                bottom_annotation = ta)
      print("HERE")
    }else{
      H=Heatmap(t(h[rownames(these_samples_metadata_rn),]),
                cluster_rows=cluster_rows,
                cluster_columns=cluster_columns,
                bottom_annotation = ta)
    }
  }

  return(H)
}


#' Assign a colour palette to metadata columns automatically and consistently
#'
#' @param metadataColumns Names of the metadata columns to assign colours for
#' @param these_samples_metadata Metadata for just the samples you need colours for
#' @param column_alias A list of column_names with values indicating what gambl colour scheme they should use (e.g. pos_neg, pathology, lymphgen)
#' @param annoAlpha Optional alpha to apply to annotation colours
#' @return Either a vector or list of colours
#' @export
#' @import dplyr ggsci
#'
#' @examples
#' all_cols=map_metadata_to_colours(legend_metadata_columns,these_meta,verbose=T)
#' all_cols=map_metadata_to_colours(c("lymphgen","pathology","genetic_subgroup"),these_samples_metadata = all_meta,column_alias=list("nothing"="FL"),as_vector = F)
map_metadata_to_colours = function(metadataColumns,
                                   these_samples_metadata,
                                   column_alias=list(),
                                   as_vector=TRUE,
                                   verbose=F,annoAlpha=1){

  clinical_colours = ggsci::get_ash("clinical")
  all_gambl_colours=get_gambl_colours()
  colours = list()
  colvec = c()
  #aliases for finding specific sets of colours
  aliases = list("COO_consensus"="coo","COO"="coo","DHITsig_consensus"="coo",
                 "pathology"="pathology","analysis_cohort"="pathology",
                 "group"="pathology",
                 "FL_group"="FL",
                 "lymphgen"="lymphgen",
                 "lymphgen_with_cnv"="lymphgen",
                 "bcl2_ba"="pos_neg",
                 "BCL2_status"="pos_neg",
                 "myc_ba"="pos_neg","bcl6_ba"="pos_neg","manta_BCL2_sv"="pos_neg",
                 "manual_BCL2_sv"="pos_neg",
                 "manta_MYC_sv"="pos_neg"
                 )
  aliases=c(aliases,column_alias)
  #print(aliases)
  for(column in metadataColumns){

    this_value = these_samples_metadata[[column]]
  options = this_value
    if(verbose){
      print(">>>>>>>")
      message("finding colour for",this_value)
      print("<<<<<<<")
    }
    if(column %in% names(aliases)){

      key = aliases[[column]]
      if(verbose){
        print(paste("using alias to look up colours for",column))
        message(paste("using",key,"for",column))
      }
      these = get_gambl_colours(classification = key)
      colours[[column]]=these
      colvec = c(colvec,these[this_value])
      if(verbose){
        message("adding:",these[this_value])
      }
    }else if(column == "sex"){
      these = get_gambl_colours("sex",alpha=annoAlpha)
      #these = these[levels(these)]
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
      colvec = c(colvec,these[this_value])
      message("adding:",these[this_value])
    }else if(sum(levels(options) %in% names(clinical_colours))==length(levels(options))){
      #we have a way to map these all to colours!
      if(verbose){
        message(paste("found colours for",column,"in clinical"))
      }
      these = clinical_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
      colvec = c(colvec,these[this_value])
    }else if(("positive" %in% options | "POS" %in% options) & length(options)<4){
      if(verbose){
        print("using pos_neg")
      }

      these = get_gambl_colours("pos_neg",alpha=annoAlpha)
      these = these[levels(options)]


      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")

      }
      colours[[column]]=these
      colvec = c(colvec,these[this_value])
    }else if("GCB" %in% options){
      these=get_gambl_colours("COO",alpha=annoAlpha)
      #names(these)=levels(options)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")

      }
      colours[[column]]=these
      colvec = c(colvec,these)
    }else if(column %in% c("pathology")){
      these=get_gambl_colours(column,alpha=annoAlpha)

      #names(these)=levels(options)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
      colvec = c(colvec,these)
    }else if(grepl("lymphgen",column,ignore.case = TRUE)){
      these=get_gambl_colours("lymphgen",alpha=annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
      colvec = c(colvec,these)
    }else if(column == "HMRN"){
      these=get_gambl_colours("hmrn",alpha=annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
      colvec = c(colvec,these)
    }else if(sum(levels(options) %in% names(all_gambl_colours))==length(levels(options))){
      if(verbose){
        message(paste("found colours for",column,"in all_gambl_colours"))
      }
      these = all_gambl_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
      colvec = c(colvec,these)
    }else if(length(levels(options))>15){

      these=rainbow(length(levels(options)),alpha=annoAlpha)
      names(these)=levels(options)

      colours[[column]]=these
      colvec = c(colvec,these)
    }else{
      these = blood_cols[sample(c(1:length(blood_cols)),size=length(levels(options)))]
      names(these)=levels(options)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
      colvec = c(colvec,these)
    }
  }
  if(as_vector){
    return(colvec)
  }
  return(colours)
}



#' Plot a sample-centric circos overview
#'
#' @param sample_id Sample ID for the sample to plot
#' @param sv_df Optional data frame of SVs (default is to use the database)
#' @param cnv_df Optional data frame of CNVs (default is to use the database)
#' @param include_sv Default TRUE
#' @param include_cnv Default TRUE
#'
#' @return Nothing
#' @export
#' @import circlize ComplexHeatmap
#'
#' @examples
#' this_samp = "13-38657_tumorB"
#' GAMBLR::plot_sample_circos(this_sample_id=this_samp,legend_metadata_columns = c("pathology","lymphgen","COO_consensus","DHITsig_consensus","bcl2_ba","myc_ba"),legend_metadata_names = c("pathology","LymphGen","COO","DHITsig","BCL2","MYC"),auto_label_sv = TRUE,chrom_list = c("chr2","chr3","chr8","chr14","chr18"))
plot_sample_circos = function(this_sample_id,sv_df,cnv_df,ssm_df,
                              include_sv=TRUE,include_ssm=FALSE,
                              legend_metadata_columns,
                              legend_metadata_names=c(),
                              include_cnv=TRUE,
                              chrom_list,label_genes,auto_label_sv=FALSE){

  add_cnv= function(cnv_df){

    bed = data.frame(cnv_df[,c("chrom","start","end","log.ratio")])
    colnames(bed)=c("chr","start","end","value1")
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = col_fun(value),
                         border = NA, posTransform = NULL, ...)
      i = getI(...)
      cell.xlim = get.cell.meta.data("cell.xlim")
    }, bg.border = NA)
  }
  if(missing(cnv_df)){
    cnv_df = get_sample_cn_segments(sample_id = this_sample_id,with_chr_prefix = TRUE)

  }
  if(missing(chrom_list)){
   chrom_list = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
  }
  if(!missing(label_genes)){
    gene_bed = grch37_gene_coordinates %>%
      dplyr::filter(gene_name %in% label_genes) %>%
      dplyr::select(chromosome,start,end,gene_name) %>%
      dplyr::mutate(chromosome = paste0("chr",chromosome))
  }
  if(missing(sv_df)){
    sv_df = get_manta_sv(with_chr_prefix = TRUE) %>%
      dplyr::filter(tumour_sample_id == this_sample_id)

  }

  sv_df = sv_df %>% dplyr::filter(CHROM_A %in% chrom_list) %>%
    dplyr::filter(CHROM_B %in% chrom_list)
  if(auto_label_sv){
    # annotate oncogene SVs and label them
    annotated_sv  = annotate_sv(sv_df,with_chr_prefix = TRUE) %>%
      dplyr::filter(!is.na(partner)) %>%
      dplyr::filter(tumour_sample_id == this_sample_id)
    these_oncogenes=unique(pull(annotated_sv,gene))
    these_partners = unique(pull(annotated_sv,partner))
    if("IGH" %in% these_partners){
      these_partners = c(these_partners,"IGHV3-62")
    }
    anno_bed1 = annotated_sv %>% dplyr::select(chrom1,start1,end1,tumour_sample_id)
    anno_bed2 = annotated_sv %>% dplyr::select(chrom2,start2,end2,tumour_sample_id)
    colnames(anno_bed1)=c("chrom","start","end","sample_id")
    colnames(anno_bed2)=c("chrom","start","end","sample_id")

    bed_mut_partner = grch37_partners %>%
      dplyr::filter(gene %in% these_partners) %>%
      mutate(chrom=paste0("chr",chrom))
    bed_mut_onco =
      grch37_oncogene %>%
      dplyr::filter(gene %in% these_oncogenes) %>%
      mutate(chrom=paste0("chr",chrom))
    bed_mut = bind_rows(bed_mut_partner,bed_mut_onco)
    print(bed_mut)
  }
  bed1 = sv_df %>% dplyr::select(CHROM_A,START_A,END_A,tumour_sample_id)
  bed2 = sv_df %>% dplyr::select(CHROM_B,START_B,END_B,tumour_sample_id)
  colnames(bed1)=c("chrom","start","end","sample_id")
  colnames(bed2)=c("chrom","start","end","sample_id")
  circos.clear()
  circos.initializeWithIdeogram(chromosome.index = chrom_list)
  add_cnv(cnv_df)
  circos.genomicLink(bed1, bed2,col = "#bdbdc1")
  if(!missing(label_genes)){
    circos.genomicLabels(gene_bed,labels.column="gene_name")
  }
  if(auto_label_sv){
    circos.genomicLink(anno_bed1, anno_bed2,col = 'red')
    circos.genomicLabels(bed_mut,labels.column="gene")
  }
  text(c(0.75,0.75),this_sample_id,cex=0.8)
  if(!missing(legend_metadata_columns)){
    samp_meta = get_gambl_metadata() %>%
      dplyr::filter(sample_id == this_sample_id)
    these_meta = samp_meta[legend_metadata_columns]
    these_cols = get_gambl_colours()
    vals = as.character(these_meta)
    names = colnames(these_meta)

    all_cols=map_metadata_to_colours(legend_metadata_columns,these_meta,verbose=T)

    cols=all_cols[vals]
    print(cols)
    if(length(legend_metadata_names) == length(legend_metadata_columns)){
      for(i in c(1:length(vals))){
        if(!legend_metadata_names[i]==""){
          vals[i]=paste(legend_metadata_names[i],vals[i])
        }
      }
    }


    lgd_discrete = Legend(labels = vals,
                          title_position = "topleft",
                          legend_gp = gpar(fill = cols))
    #draw(lgd_discrete, x = unit(4, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
    draw(lgd_discrete, x = unit(3, "mm"), y = unit(1, "npc") - unit(5,"mm"),
         just = c("left", "top"))
  }

  # continuous
  lgd_cnv = Legend(at = c(-2, -1, 0, 1, 2),
                     col_fun=circlize::colorRamp2(c(-1, 0, 1),
                                    c("blue", "white", "red")),
                     title_position = "topleft", title = "log\nratio")

  draw(lgd_cnv, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"),
       just = c("right", "bottom"))
}

#' Make an oncoplot that is pretty using ComplexHeatmap. The metadata is expected to follow the structure and column naming used in GAMBL.
#' If you provide your own non-GAMBL samples and metadata, you must include at least the following columns with these names.
#' The first one should match the Tumor_Sample_Barcode in the MAF object or onco_matrix you provide.
#' sample_id, pathology
#'
#' @param maftools_obj A maftools object containing the mutations you want to plot
#' @param onco_matrix_path Provide a path to an onco_matrix file instead of a MAF object if the former is unavailable (this limits functionality a bit)
#' @param genes An optional list of genes to restrict your plot to
#' @param keepGeneOrder Set to TRUE if you want to preserve the gene order specified
#' @param keepSampleOrder Set to TRUE if you want to preserve the sample order specified
#' @param these_samples_metadata Data frame containing metadata for your samples
#' @param metadataColumns A vector containing the categorical column names you want to plot below
#' @param numericMetadataColumns A vector containing the numeric columns you want to plot below
#' @param numericMetadataMax A numeric vector of cutoffs to apply to numeric columns above
#' @param sortByColumns A vector containing the column names you want to sort columns (patients) on
#' @param removeNonMutated Set to TRUE to drop unmutated cases
#' @param minMutationPercent Only genes mutated in more than minMutationPercent % patients will be included
#' @param fontSizeGene Font size for gene labels (default 6)
#' @param annoAlpha Optional alpha to apply to annotation colours
#' @param mutAlpha Optional alpha to apply to mutation colours
#' @param recycleOncomatrix Set to TRUE most of the time to reuse the oncomatrix saved by maftools
#' @param box_col Colour of boxes for outlining mutations (can be problematic with larger oncoprints)
#' @param legend_row Fiddle with these to widen or narrow your legend
#' @param legend_col Fiddle with these to widen or narrow your legend
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette
#'
#' @return
#' @export
#' @import ComplexHeatmap grid
#'
#' @examples
#' prettyOncoplot(maftools_obj = maf_obj,genes = bl_genes,
#' these_samples_metadata = extra_meta,
#' metadataColumns = c("pathology","COO_consensus",
#'                     "cluster_name", "lymphgen","EBV_status_inf",
#'                     "manta_BCL6_sv"),
#' hide_annotations = c(some_ashm,"lymphgen","COO_consensus",
#'                      "pathology","manta_BCL6_sv"),
#' expressionColumns = c("IRF4",some_ashm),
#' sortByColumns = c("IRF4"),
#' keepGeneOrder = FALSE,splitGeneGroups = groups,
#' splitColumnName = "cluster_name",
#' metadataBarHeight = 2.5,metadataBarFontsize = 6,fontSizeGene = 8,
#' recycleOncomatrix = TRUE,removeNonMutated = FALSE)
prettyOncoplot = function(maftools_obj,
                          onco_matrix_path,
                          genes,
                          include_noncoding=list("NFKBIZ"=c("3'UTR"),"HNRNPH1"="Splice_Region"),
                          keepGeneOrder=TRUE,
                          keepSampleOrder=TRUE,
                          highlightHotspots=FALSE,
                          these_samples_metadata,
                          metadataColumns,
                          numericMetadataColumns,
                          expressionColumns=c(),
                          numericMetadataMax,sortByColumns,
                          removeNonMutated=FALSE,
                          minMutationPercent,
                          fontSizeGene=6,
                          annoAlpha=1,mutAlpha=1,
                          recycleOncomatrix=FALSE,
                          box_col=NA,
                          metadataBarHeight=1.5,
                          metadataBarFontsize=5,
                          hideTopBarplot=TRUE,
                          hideSideBarplot=FALSE,
                          splitColumnName,
                          splitGeneGroups,
                          legend_row=3,
                          legend_col=3,
                          showTumorSampleBarcode=FALSE,
                          groupNames,verbose=FALSE,
                          hide_annotations,
                          custom_colours = NULL,
                          legend_direction="horizontal",
                          legend_position="bottom",
                          annotation_row=2,
                          annotation_col=1,
                          legendFontSize=10
                          ){

  patients = pull(these_samples_metadata,sample_id)
  #ensure patients not in metadata get dropped up-front to ensure mutation frequencies are accurate
  if(!recycleOncomatrix & missing(onco_matrix_path)){
    onco_matrix_path="onco_matrix.txt"
  #order the data frame the way you want the patients shown
    maf_patients =unique(as.character(maftools_obj@data$Tumor_Sample_Barcode))
    if(any(!maf_patients %in% patients)){
      extra = maf_patients[which(maf_patients %in% patients)]
      patients = maf_patients[which(maf_patients %in% patients)]
      n_drop=length(extra)
      message(paste(n_drop,"patients are not in your metadata, will drop them from the data before displaying"))
      maftools_obj = subsetMaf(maf = maftools_obj,tsb=patients)

    }
    if(missing(genes)){
      #check that our MAFtools object only contains samples in the supplied metadata

      genes = maftools::getGeneSummary(x = maftools_obj)[order(MutatedSamples, decreasing = TRUE)][,.(Hugo_Symbol, MutatedSamples)]
      colnames(genes)[2] = "mutload"
      totSamps = as.numeric(maftools_obj@summary[3,summary])
      genes[,fractMutated := mutload/totSamps]

      genes = genes[fractMutated*100 >= minMutationPercent, Hugo_Symbol]

      lg = length(genes)
      message(paste("creating oncomatrix with",lg,"genes"))
      om = maftools:::createOncoMatrix(m = maftools_obj,
                                     g = genes,
                                     add_missing = TRUE)
      mat_origin = om$oncoMatrix
      tsbs = levels(maftools:::getSampleSummary(x = maftools_obj)[,Tumor_Sample_Barcode])
      print(paste("numcases:",length(tsbs)))
      if(!removeNonMutated){
        tsb.include = matrix(data = 0, nrow = nrow(mat_origin),
                             ncol = length(tsbs[!tsbs %in% colnames(mat_origin)]))
        colnames(tsb.include) = tsbs[!tsbs %in% colnames(mat_origin)]
        rownames(tsb.include) = rownames(mat_origin)

        mat_origin = cbind(mat_origin, tsb.include)
      }
      write.table(mat_origin,file=onco_matrix_path,quote=F,sep="\t")

    }else{
      om = maftools:::createOncoMatrix(m = maftools_obj,
                                       g = genes,
                                       add_missing = TRUE)
      mat_origin = om$oncoMatrix
      tsbs = levels(maftools:::getSampleSummary(x = maftools_obj)[,Tumor_Sample_Barcode])
      print(paste("numcases:",length(tsbs)))
      if(!removeNonMutated){
        tsb.include = matrix(data = 0, nrow = nrow(mat_origin),
                             ncol = length(tsbs[!tsbs %in% colnames(mat_origin)]))
        colnames(tsb.include) = tsbs[!tsbs %in% colnames(mat_origin)]
        rownames(tsb.include) = rownames(mat_origin)

        mat_origin = cbind(mat_origin, tsb.include)
      }
      write.table(mat_origin,file=onco_matrix_path,quote=F,sep="\t")
      #maftools::oncoplot(maftools_obj,genes=genes,writeMatrix = T,removeNonMutated = removeNonMutated)
    }
  }
  if(missing(onco_matrix_path)){
    onco_matrix_path="onco_matrix.txt"
  }
  if(!missing(numericMetadataColumns)){
    message(paste0("The column(s) ", numericMetadataColumns, " specified both in metadata and numeric metadata. Plotting as numeric values..."))
    metadataColumns = metadataColumns[!metadataColumns %in% numericMetadataColumns]
  }
  patients = pull(these_samples_metadata,sample_id)
  #because the way MAFtools writes this file out is the absolute worst for compatability
  old_style_mat = read.table(onco_matrix_path,sep="\t",stringsAsFactors = FALSE)
  mat=read.table(onco_matrix_path,sep="\t",header=TRUE,check.names = FALSE,row.names=1,fill=TRUE,stringsAsFactors = F,na.strings = c("NA",""))
  colnames(old_style_mat) = colnames(mat)
  mat=old_style_mat
  mat[mat==0]=""
  #add the noncoding mutations to this if requested (only for genes and types specified)
  if(length(include_noncoding)>0){
    all_genes_df = data.frame(Hugo_Symbol=rownames(mat))
    all_samples_df = data.frame(Tumor_Sample_Barcode=colnames(mat))
    for(gene in names(include_noncoding)){
      for(this_vc in unname(include_noncoding[[gene]])){
        message(paste(gene,"and",this_vc))
        these_samples = dplyr::filter(maftools_obj@maf.silent,Hugo_Symbol == gene &
                                      Variant_Classification == this_vc) %>%
        dplyr::select(Tumor_Sample_Barcode,Variant_Classification) %>%
        unique() %>% pull(Tumor_Sample_Barcode)
        for(samp in these_samples){
          if(samp %in% colnames(mat)){
            if(mat[gene,samp]==""){
              mat[gene,samp] = this_vc
            }else{
              mat[gene,samp] = paste0(this_vc,";",mat[gene,samp])
            }
          }
        }
      }
    }
  }
  #annotate hot spots if necessary
  if(missing(metadataColumns)){
    message("you should name at least one metadata column to show as an annotation. Defaulting to pathology")
    metadataColumns=c("pathology")
  }

  if(missing(genes)){
    genes = rownames(mat)
  }

  col=get_gambl_colours("mutation",alpha=mutAlpha)
  mat[mat==0]=""


  patients_kept = patients[which(patients %in% colnames(mat))]
  patients_dropped = patients[which(!patients %in% colnames(mat))]
  if(verbose){
    message("====DROPPED=====")
    message(patients_dropped)
  }

  genes_kept = genes[which(genes %in% rownames(mat))]
  if(!missing(minMutationPercent)){
    if(! onco_matrix_path == "onco_matrix.txt"){

      warning("mintMutationPercent option is not available when you provide your own oncomatrix. Feel free to implement this if you need it")
      return()
    }

    mutation_counts = maftools_obj@gene.summary %>%
      mutate(fake_column=1) %>%
      complete(., tidyr::expand(., crossing(fake_column), Hugo_Symbol = genes)) %>%
      select(-fake_column) %>%
      replace(is.na(.), 0) %>%
      dplyr::select(Hugo_Symbol,MutatedSamples) %>%
      as.data.frame()
    numpat=length(patients)
    mutation_counts = mutate(mutation_counts,percent_mutated=100 * MutatedSamples/numpat)
    genes_keep = mutation_counts %>%
      dplyr::filter(percent_mutated>=minMutationPercent) %>%
      pull(Hugo_Symbol)
    genes_kept = genes[genes %in% genes_keep]
  }

  mat = mat[,patients_kept]
  mat = mat[which(rownames(mat) %in% genes_kept),]
  spacing = 0
  height_scaling = 1

  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#e6e6e6", col = box_col))
    },
    RNA = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#F2ED36", col = box_col))
    },
    `3'UTR` = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#F2ED36", col = box_col))
    },
    Intron = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), 0.75* h*height_scaling,
                gp = gpar(fill = col["Nonsense_Mutation"], col = box_col))
    },
    # big blue
    Nonsense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#D8A7CA", col = box_col))
    },
    Splice_Site = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Splice_Site"], col = box_col))
    },
    Splice_Region = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Splice_Region"], col = box_col))
    },
    Nonstop_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Nonstop_Mutation"], col = box_col))
    },
    Translation_Start_Site = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Translation_Start_Site"], col = box_col))
    },
    In_Frame_Ins = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["In_Frame_Ins"], col = box_col))
    },
    In_Frame_Del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["In_Frame_Del"], col = box_col))
    },
    #all frame shifts will be the same colour, magenta
    Frame_Shift_Del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Frame_Shift_Del"], col = box_col))
    },
    Frame_Shift_Ins = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Frame_Shift_Ins"], col = box_col))
    },
    # big red
    Multi_Hit = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Multi_Hit"], col = box_col))
    },
    # small green
    Missense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Missense_Mutation"], col = box_col))
    },
    hot_spot = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), (height_scaling/5)*h,
                gp = gpar(fill = "white", col = box_col))
    }
  )
  #automagically assign colours for other metadata columns.
  #TO DO: convert the loop below into a "map_metadata_to_colours" function
  blood_cols=get_gambl_colours("blood",alpha=annoAlpha)

  colours = list()
  clinical_colours = ggsci::get_ash("clinical")
  all_gambl_colours=get_gambl_colours()

  for(column in metadataColumns){

    these_samples_metadata[[column]] = factor(these_samples_metadata[[column]], levels=unique(these_samples_metadata[[column]]))
    options = these_samples_metadata %>% arrange(column) %>% dplyr::filter(!is.na(column)) %>% pull(column) %>% unique()
    options=options[!is.na(options)]
    if(verbose){
      print(">>>>>>>")
      print(levels(options))
      print("<<<<<<<")
    }
    if(column == "sex"){
      these = get_gambl_colours("sex",alpha=annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
    }else if(sum(levels(options) %in% names(clinical_colours))==length(levels(options))){
      #we have a way to map these all to colours!
      if(verbose){
        message(paste("found colours for",column, "here"))
      }
      these = clinical_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
    }else if(("positive" %in% options | "POS" %in% options | "yes" %in% options) & length(options)<4){
      if(verbose){
        print("using pos_neg")
      }

      these = get_gambl_colours("pos_neg",alpha=annoAlpha)
      these = these[levels(options)]

      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")

      }
      colours[[column]]=these
    }else if("GCB" %in% options){
      these=get_gambl_colours("COO",alpha=annoAlpha)
      #names(these)=levels(options)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")

      }
      colours[[column]]=these
    }else if(column %in% c("pathology")){
      these=get_gambl_colours(column,alpha=annoAlpha)

      #names(these)=levels(options)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
    }else if(grepl("lymphgen",column,ignore.case = TRUE)){
      these=get_gambl_colours("lymphgen",alpha=annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
    }else if(column == "HMRN"){
      these=get_gambl_colours("hmrn",alpha=annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these

    }else if(sum(levels(options) %in% names(all_gambl_colours))==length(levels(options))){
      if(verbose){
        message(paste("found colours for",column,"in all_gambl_colours"))
      }
      these = all_gambl_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
    }else if(length(levels(options))>15){

      these=rainbow(length(levels(options)),alpha=annoAlpha)
      names(these)=levels(options)

        colours[[column]]=these

    }else{
      these = blood_cols[sample(c(1:length(blood_cols)),size=length(levels(options)))]
      names(these)=levels(options)
      if(!"NA" %in% names(these)){
        these= c(these,"NA"="white")
      }
      colours[[column]]=these
    }
  }

  if (! is.null(custom_colours)){
    colours=custom_colours
  }

  if(highlightHotspots){
    hot_samples = dplyr::filter(maftools_obj@data,hot_spot==TRUE & Hugo_Symbol %in% genes) %>%
    dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode) %>% mutate(mutated="hot_spot") %>% unique()
    all_genes_df = data.frame(Hugo_Symbol=rownames(mat))
    all_samples_df = data.frame(Tumor_Sample_Barcode=colnames(mat))
    hs = left_join(all_samples_df,hot_samples)
    hot_mat = hs %>%
        pivot_wider(names_from = "Tumor_Sample_Barcode",values_from="mutated") %>%
    left_join(all_genes_df,.) %>%
    column_to_rownames("Hugo_Symbol") %>% as.matrix()
    # annotate hotspots in matrix
    for (i in colnames(mat)){
      mat[genes,i][!is.na(hot_mat[genes,i])] <- paste0(mat[genes,i][!is.na(hot_mat[genes,i])], ";",  hot_mat[genes,i][!is.na(hot_mat[genes,i])])
    }
    colours[["hot_spots"]]= c("hot_spot"="magenta")
  }





  if(verbose){
    print(colours) #eventually get rid of this once the bugs are gone
  }
  if(!missing(numericMetadataColumns)){
    metadata_df = dplyr::filter(these_samples_metadata, sample_id %in% patients_kept) %>%
      column_to_rownames("sample_id") %>%
      dplyr::select(all_of(c(metadataColumns,numericMetadataColumns,expressionColumns)))
    if(!missing(numericMetadataMax)){
      max_list = setNames(numericMetadataMax,numericMetadataColumns)

      metadata_df = metadata_df %>%
        mutate(across(names(max_list), ~ ifelse(.x > max_list[[cur_column()]], max_list[[cur_column()]], .x)))
    }

  }else{
    metadata_df = dplyr::filter(these_samples_metadata, sample_id %in% patients_kept) %>%
    column_to_rownames("sample_id") %>% dplyr::select(all_of(c(metadataColumns,expressionColumns)))
  }
  metadata_df = metadata_df %>%
    mutate(across(all_of(expressionColumns), ~ trim_scale_expression(.x)))

  if(!missing(sortByColumns)){
    metadata_df = arrange(metadata_df,across(sortByColumns))
    patients_kept = rownames(metadata_df)
  }
  if(verbose){
    print(genes_kept)
  }
  col_fun=circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  for(exp in expressionColumns){
    colours[[exp]] = col_fun
  }

  if(missing(splitColumnName)){
    column_split=rep("",length(patients_kept))
  }else{
    column_split = factor(metadata_df[patients_kept,splitColumnName])
  }
  if(missing(splitGeneGroups)){
    row_split=rep("",length(genes))
  }else{
    row_split=factor(splitGeneGroups[genes],levels=unique(splitGeneGroups[genes]))
  }
  if(!missing(groupNames)){
    column_title=groupNames
  }else{
    column_title=NULL
  }
  if(keepGeneOrder){
    gene_order=genes
  }else{
    gene_order=NULL
  }
  if(missing(hide_annotations)){
    show_legend = rep(TRUE,length(colnames(metadata_df)))
  }else{
    show_legend = rep(TRUE,length(colnames(metadata_df)))
    names(show_legend) = colnames(metadata_df)
    show_legend[hide_annotations]=FALSE
  }
  if(missing(sortByColumns)){
    column_order = NULL
  }else{
    column_order = patients_kept
  }
  heatmap_legend_param = list(title = "Alterations",
                         at = c("RNA", "3'UTR" , "Nonsense_Mutation", "Splice_Site","Splice_Region", "Nonstop_Mutation", "Translation_Start_Site",
                         "In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Frame_Shift_Del", "Multi_Hit", "Missense_Mutation", "hot_spot"), 
                         labels = c("RNA", "3'UTR", "Nonsense Mutation", "Splice Site","Splice Region", "Nonstop Mutation", "Translation Start Site",
                         "In Frame Insertion", "In Frame Deletion", "Frame Shift Insertion", "Frame Shift Deletion",
                         "Multi Hit", "Missense Mutation", "Hotspot"),
                         nrow=annotation_row, ncol=annotation_col,
                         legend_direction = legend_direction,
                         labels_gp = gpar(fontsize = legendFontSize))
  if(hideTopBarplot){
    top_annotation = NULL
  }else{
    top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot())
  }

  ch = ComplexHeatmap::oncoPrint(mat[intersect(genes, genes_kept),patients_kept],
                                   alter_fun = alter_fun,
                                   top_annotation=top_annotation,
                                   right_annotation=NULL,
                                   col = col,
                                   row_order=gene_order,
                                   column_order = column_order,
                                   column_labels = NULL,
                                   show_column_names = showTumorSampleBarcode,
                                   column_split=column_split,
                                   column_title=column_title,
                                   row_title=NULL,
                                   row_split=row_split[intersect(genes, genes_kept)],
                                   heatmap_legend_param = heatmap_legend_param,
                                   row_names_gp = gpar(fontsize = fontSizeGene),
                                   pct_gp = gpar(fontsize = fontSizeGene),
                    bottom_annotation =
                    ComplexHeatmap::HeatmapAnnotation(df=metadata_df,
                                                      show_legend=show_legend,
                                                      col=colours,
                                                      simple_anno_size = unit(metadataBarHeight, "mm"),
                                                      gap = unit(0.25*metadataBarHeight, "mm"),
                                                      annotation_name_gp=gpar(fontsize=metadataBarFontsize),
                                                      annotation_legend_param =
                                                        list(nrow=legend_row,
                                                            col_fun=col_fun,
                                                            ncol=legend_col,
                                                            direction=legend_direction,
                                                            labels_gp = gpar(fontsize = legendFontSize))))
    draw(ch, heatmap_legend_side = legend_position, annotation_legend_side = legend_position)
}

#' Generate a colourful multi-panel overview of hypermutation in regions of interest across many samples
#'
#' @param regions_bed Bed file with chromosome coordinates, should contain columns chr, start, end, name (with these exact names)
#' @param regions_to_display Optional vector of names from default regions_bed to use
#' @param metadata A metadata file already subsetted and arranged on the order you want the samples vertically displayed
#' @param classification_column optional. Override default column for assigning the labels used for colouring in the figure.
#'
#' @return nothing
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
ashm_multi_rainbow_plot = function(regions_bed,regions_to_display,
                                   exclude_classifications,metadata,custom_colours,
                                   classification_column="lymphgen",maf_data){
  table_name = config::get("results_tables")$ssm
  db=config::get("database_name")
  #get the mutations for each region and combine
  #regions_bed should contain chr, start, end, name (with these exact names)
  if(missing(metadata)){
    metadata = get_gambl_metadata()
    meta_arranged = arrange(metadata,pathology_rank,lymphgen)
  }else{
    meta_arranged = metadata #assume the user already arranged it the way they wanted
  }
  if(!missing(exclude_classifications)){
    meta_arranged = dplyr::filter(meta_arranged,!get(classification_column) %in% exclude_classifications)
  }
  if(missing(regions_bed)){
    regions_bed= grch37_ashm_regions
    regions_bed = mutate(regions_bed,regions=paste0(chr_name,":",hg19_start,"-",hg19_end))
    regions_bed = mutate(regions_bed,name=paste0(gene,"-",region))
  }else{
    regions_bed = mutate(regions_bed,regions=paste0(chr,":",start,"-",end))
    #if name column is missing, add it
    if(!"name" %in% colnames(regions_bed))
    {
      regions_bed$name =regions_bed$regions
    }
  }
  print(regions_bed)
  names=pull(regions_bed,name)
  names = c(names,"NFKBIZ-UTR","MAF","PAX5","WHSC1","CCND1",
            "FOXP1-TSS1","FOXP1-TSS2","FOXP1-TSS3","FOXP1-TSS4","FOXP1-TSS5",
            "BCL6","IGH","IGL","IGK","PVT1","BCL2") #add some additional regions of interest
  regions = pull(regions_bed,regions)
  regions = c(regions,"chr3:101578214-101578365","chr16:79627745-79634622","chr9:36898851-37448583","chr4:1867076-1977887","chr11:69451233-69460334","chr3:71623481-71641671","chr3:71532613-71559445","chr3:71343345-71363145","chr3:71167050-71193679","chr3:71105715-71118362",
              "chr3:187406804-188522799","chr14:106144562-106344765","chr22:23217074-23250428","chr2:89073691-89320640",
              "chr8:128774985-128876311","chr18:60982124-60990180")
  regions_bed = data.frame(regions=regions,names=names)
  regions_bed = dplyr::filter(regions_bed,names %in% regions_to_display)
  regions = pull(regions_bed,regions)
  names=pull(regions_bed,names)

  if(missing(maf_data)){
    region_mafs = lapply(regions,function(x){get_ssm_by_region(region=x,streamlined = TRUE)})
  }else{
    region_mafs = lapply(regions,function(x){get_ssm_by_region(region=x,streamlined=TRUE,maf_data=maf_data)})
  }

  tibbled_data = tibble(region_mafs, region_name = names)
  unnested_df = tibbled_data %>% unnest_longer(region_mafs)
  unlisted_df = mutate(unnested_df,start=region_mafs$Start_Position,sample_id=region_mafs$Tumor_Sample_Barcode) %>%
    dplyr::select(start,sample_id,region_name)



  meta_arranged$classification = factor(meta_arranged[,classification_column],levels=unique(meta_arranged[,classification_column]))
  muts_anno = left_join(unlisted_df,meta_arranged)
  muts_first =  dplyr::select(muts_anno,start,region_name) %>% group_by(region_name) %>% arrange(start) %>% dplyr::filter(row_number()==1)
  eg = expand_grid(start=pull(muts_first,start),sample_id=pull(meta_arranged,sample_id))
  eg = left_join(eg,muts_first)

  #concatenate expanded frame of points with original mutation data
  real_and_fake = bind_rows(unlisted_df,eg)
  muts_anno = left_join(real_and_fake,meta_arranged) #%>% filter(!is.na(get(classification_column)))

  muts_anno$sample_id= factor(muts_anno$sample_id,levels=meta_arranged$sample_id)

  if(!missing(regions_to_display)){
    muts_anno = dplyr::filter(muts_anno,region_name %in% regions_to_display)
  }
  #make the plot
  if(missing(custom_colours)){
    p = muts_anno %>%
      ggplot() + geom_point(aes(x=start,y=sample_id,colour=classification),alpha=0.4,size=0.6) +
      theme(axis.text.y=element_blank()) +
      facet_wrap(~region_name,scales="free_x") + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))

    print(p)
  }else{
    #testing manual colouring
    p = muts_anno %>%
      ggplot() + geom_point(aes(x=start,y=sample_id,colour=classification),alpha=0.4,size=0.6) +
      theme(axis.text.y=element_blank()) + scale_colour_manual(values=custom_colours) +
      facet_wrap(~region_name,scales="free_x") + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
    print(p)
  }
}

#' Create a genome-wide copy number plot for one sample and (optionally) display mutation VAF
#'
#' @param this_sample The sample_id for the sample to plot
#' @param just_segments Specify whether only the segments will be plotted (instead of mutation VAF)
#' @param genes_to_label optional. Provide a list of genes to label (if mutated). Can only be used with coding_only (see below)
#' @param coding_only optional. Set to TRUE to restrict to plotting only coding mutations
#' @param from_flatfile If set to true the function will use flatfiles instead of the database
#'
#' @return nothing
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
copy_number_vaf_plot = function(this_sample,just_segments=FALSE,coding_only=FALSE,one_chrom,
                                genes_to_label,from_flatfile=FALSE,use_augmented_maf=FALSE,add_chr_prefix = FALSE){
  chrom_order=factor(c(1:22,"X"))
  if(add_chr_prefix){
    chrom_order = c(1:22,"X")
    chrom_order = factor(unlist(lapply(chrom_order,function(x){paste0("chr",x)})))
  }
  cn_colours = get_gambl_colours(classification = "copy_number")
  maf_and_seg = assign_cn_to_ssm(this_sample=this_sample,coding_only=coding_only,from_flatfile=from_flatfile,use_augmented_maf=use_augmented_maf)
  vaf_cn_maf = maf_and_seg[["maf"]]
  vaf_cn_maf = mutate(vaf_cn_maf,CN=case_when(LOH == "1" & CN == 2 ~ "nLOH",
                                              TRUE ~ as.character(CN)))
  if(!missing(one_chrom)){
    vaf_cn_maf = dplyr::filter(vaf_cn_maf,Chromosome == one_chrom)
  }
  #vaf_cn_maf = mutate(vaf_cn_maf,CN=as.character(CN))
  #use_neut_loh=TRUE

  if(just_segments){
    #I realized this is ugly
    cn_seg = maf_and_seg[["seg"]]
    cn_seg = mutate(cn_seg,CN_segment = as.numeric(CN),CN = as.character(CN))
    print(head(cn_seg))
    if(!missing(one_chrom)){
      cn_seg = dplyr::filter(cn_seg,Chromosome %in% one_chrom)
    }

    ggplot(cn_seg) +
      geom_segment(data=cn_seg,aes(x=Start_Position,xend=End_Position,y=CN_segment,yend=CN_segment,colour=CN)) +
      facet_wrap(~Chromosome,scales="free_x") +
      scale_colour_manual(values = cn_colours) +
      theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
  }else{
    if(coding_only){
      if(missing(genes_to_label)){
        p = mutate(vaf_cn_maf,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
          geom_point(aes(x=Start_Position,y=vaf,colour=CN),alpha=0.6,size=2) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") +
          theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
        p + ggtitle(this_sample)
      }else{
        #label any mutations that intersect with our gene list
        plot_genes = vaf_cn_maf %>% dplyr::filter(Hugo_Symbol %in% my_genes)

        p = mutate(vaf_cn_maf,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
          geom_point(aes(x=Start_Position,y=vaf,colour=CN),size=2) +
          geom_text(data=plot_genes,aes(x=Start_Position,y=0.8,label=Hugo_Symbol),size=3,angle=90) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") + ylim(c(0,1)) +
          theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
        p + ggtitle(this_sample)
      }
    }else{

      p = mutate(vaf_cn_maf,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
        geom_point(aes(x=Start_Position,y=vaf,colour=CN),alpha=0.6,size=0.2) +
        scale_colour_manual(values = cn_colours) +
        facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") +
        theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
      p + ggtitle(this_sample)
    }
  }
}


#TODO migrate viz/plotting functions that don't directly rely on the database to a separate file
#' Make a rainbow plot of all mutations in a region, ordered and coloured by metadata
#'
#' @param mutations_maf A data frame containing mutations (MAF format) within a region of interest (i.e. use the get_ssm_by_region)
#' @param metadata should be a data frame with sample_id as a column that should match Tumor_Sample_Barcode in the database
#' @param classification_column The name of the metadata column to use for ordering and colouring samples
#' @param bed Optional data frame specifying the regions to annotate (required columns: start, end, name)
#'
#' @return ggplot2 object
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' # basic usage
#' ashm_rainbow_plot(mutations_maf=get_ssm_by_region(region="chr8:128806578-128806992"),metadata=get_gambl_metadata())
#' # advanced usages
#' mybed = data.frame(start=c(128806578,128805652,128748315), end=c(128806992,128809822,128748880), name=c("TSS","enhancer","MYC-e1"))
#' ashm_rainbow_plot(mutations_maf=my_mutations,metadata=my_metadata,bed=mybed)
ashm_rainbow_plot = function(mutations_maf,
                             metadata,exclude_classifications,
                             drop_unmutated=FALSE,
                             classification_column,
                             bed,region,custom_colours,hide_ids=TRUE){
  table_name = config::get("results_tables")$ssm
  db=config::get("database_name")
  if(!missing(region)){

    region = gsub(",","",region)
    split_chunks = unlist(strsplit(region,":"))
    chromosome = split_chunks[1]

    startend = unlist(strsplit(split_chunks[2],"-"))
    qstart=as.numeric(startend[1])
    qend=as.numeric(startend[2])

    if(missing(mutations_maf)){
      mutations_maf = get_ssm_by_region(region=region,streamlined = TRUE)
    }else{
      #ensure it only contains mutations in the region specified
      mutations_maf = get_ssm_by_region(region=region,streamlined = TRUE,maf_data = mutations_maf)
    }

  }
  if(!missing(classification_column)){
    meta_arranged = arrange(metadata,pathology_rank,lymphgen) #fix this to use the actual column
    if(!missing(exclude_classifications)){
      meta_arranged = dplyr::filter(meta_arranged,!get(classification_column) %in% exclude_classifications)
    }
  }else{
    classification_column = "lymphgen"
    meta_arranged = metadata
  }


  mutation_positions = mutations_maf %>%
    dplyr::select(Tumor_Sample_Barcode,Start_Position) %>% as.data.frame()
  mutated_cases = pull(mutation_positions,Tumor_Sample_Barcode) %>% unique()
  if(drop_unmutated){
    meta_arranged = meta_arranged %>% dplyr::filter(sample_id %in% mutated_cases)
  }
  #add a fake mutation at the start position for each sample to ensure every sample shows up
  fake_mutations = data.frame(Tumor_Sample_Barcode=pull(metadata,sample_id),Start_Position=qstart-1000)
  mutation_positions = rbind(mutation_positions,fake_mutations)

  meta_arranged$classification = factor(meta_arranged[,classification_column],levels=unique(meta_arranged[,classification_column]))

  muts_anno = dplyr::left_join(mutation_positions,meta_arranged,by=c("Tumor_Sample_Barcode" = "sample_id"))


  muts_anno$sample_id = factor(muts_anno$Tumor_Sample_Barcode,levels=unique(meta_arranged$sample_id))
  #

  if(missing(custom_colours)){
    p = ggplot(muts_anno) +
      geom_point(aes(x=Start_Position,y=sample_id,colour=classification),alpha=0.4)
  }else{
    p = ggplot(muts_anno) + geom_point(aes(x=Start_Position,y=sample_id,colour=classification),alpha=0.4) +
      scale_colour_manual(values=custom_colours)
  }
  if(missing(bed)){
    p +
      guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
  }else{
    bed = bed %>% mutate(size = end - start) %>% mutate(midpoint = start + size/2)
    height = length(unique(meta_arranged$sample_id)) + 8
    p = p + geom_rect(data=bed, aes(xmin = start, xmax = end, ymin = 0, ymax = height+5),alpha=0.1) +
      geom_text(data=bed,aes(x = midpoint , y= height ,label=name),size = 2.5,angle=90) +
      guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))

  }
  if(hide_ids){
    p= p + theme(axis.text.y=element_blank())
  }else{
    p = p + theme(axis.text.y = element_text(size = 5))
  }
  return(p)
}

#' This function doesn't do anything yet
#'
#' @param mafs TODO
#' @param sample_id TODO
#' @param genes TODO
#' @param show_noncoding TODO
#' @param detail TODO
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
      coding.maf = dplyr::filter(all.maf,!Variant_Classification %in% c("Silent","RNA","IGR","Intron","5'Flank","3'Flank","5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = dplyr::filter(coding.maf,(Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }
    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = dplyr::filter(coding.maf,coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = dplyr::filter(coding.maf,coord %in% B.zero.coords)
    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords),"category"]="not-A"
    coding.maf[which(coord %in% B.zero.coords),"category"]="not-B"



    #actually this is changed in eitehr direction, not just gained
    just_gained_lg_all = dplyr::filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" )
    #just_gained_lg = filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point !=2) %>%
    #  mutate(time_point = time_point +0.4)

    just_gained_lg = dplyr::filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point == 2 & VAF > 0) %>%
      mutate(time_point = time_point +0.4)
    print(just_gained_lg)
    just_trunk = dplyr::filter(coding.maf, Hugo_Symbol %in% genes & category == "trunk" & time_point ==1) %>%
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
      coding.maf = dplyr::filter(all.maf,!Variant_Classification %in% c("Silent","RNA","IGR","Intron","5'Flank","3'Flank","5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = dplyr::filter(coding.maf,(Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }

    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = dplyr::filter(coding.maf,coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = dplyr::filter(coding.maf,coord %in% B.zero.coords)

    C.rows = which(coding.maf$time_point==3)
    C.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==3 & VAF == 0), "coord"]))
    C.zero = dplyr::filter(coding.maf,coord %in% C.zero.coords)

    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords),"category"]="not-A"
    coding.maf[which(coord %in% B.zero.coords),"category"]="not-B"
    coding.maf[which(coord %in% C.zero.coords),"category"]="not-C"

    just_gained_lg_all = dplyr::filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" )
    just_gained_lg = dplyr::filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" & time_point ==3) %>%
      mutate(time_point = time_point +0.4)

    just_trunk = dplyr::filter(coding.maf, Hugo_Symbol %in% lg & category == "trunk" & time_point ==1) %>%
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


#' Use GISTIC2.0 scores output to reproduce maftools::chromoplot with more flexibility
#'
#' @param scores output file scores.gistic from the run of GISTIC2.0
#' @param genes_to_label optional. Provide a data frame of genes to label (if mutated). The first 3 columns must contain chromosome, start, and end coordinates. Another required column must contain gene names and be named `gene`. All other columns are ignored. If no data frame provided, oncogenes from GAMBLR packages are used by default to annotate on the plot.
#' @param cutoff optional. Used to determine which regions to color as aberrant. Must be float in the range [0-1]. The higher the number, the less regions will be considered as aberrant. The default is 0.5.
#' @param adjust_amps optional. The value of G-score for highest amplification peak will be multiplied by this value to determine how far up the gene label will be displayed. Default 0.5.
#' @param adjust_dels optional. The value of G-score for highest deletion peak will be multiplied by this value to determine how far down the gene label will be displayed. Default 2.75.
#' @param label_size optional. The font size for the gene label to be displayed. Default 3.
#' @param force_pull optional. How strong the gene name label will be pulled towards a data point. Default 0 (no pulling).
#' @param segment.curvature optional. Indicates whether arrow to the data point should be curved. Accepts numeric value, where negative is for left-hand and positive for right-hand curves, and 0 for straight lines. Default 0.25
#' @param segment.ncp optional. Indicates number of control points to make a smoother curve. Higher value allows for more flexibility for the curve. Default 4
#' @param segment.angle optional. Numeric value in the range 0-180, where less than 90 skews control points of the arrow from label to data point toward the start point. Default 25
#'
#'
#' @return nothing
#' @export
#' @import tidyverse ggrepel
#'
#' @examples
#' # basic usage
#' prettyChromoplot("path_to_gistic_results/scores.gistic")
#' # advanced usages
#' prettyChromoplot("path_to_gistic_results/scores.gistic", genes_to_label="path_to_gene_coordinates_table.tsv", cutoff=0.75) +
#' ...(any ggplot options to customize plot appearance)
prettyChromoplot = function(scores,
                            genes_to_label,
                            cutoff=0.5,
                            adjust_amps=0.5,
                            adjust_dels=2.75,
                            label_size=3,
                            force_pull = 0,
                            segment.curvature = 0.25,
                            segment.ncp = 4,
                            segment.angle = 25){
  # read GISTIC scores file, convert G-score to be negative for deletions, and relocate chromosome, start, and end columns to be the first three
  scores <- data.table::fread(scores) %>%
    dplyr::mutate(`G-score`= ifelse(Type=="Amp",  `G-score`, -1*`G-score`)) %>%
    dplyr::relocate(Type, .after=frequency)
  # annotate each region with direction of changes - used for coloring
  scores$fill <- ifelse(scores$Type == "Amp" & scores$`-log10(q-value)` > cutoff, "up",
                        ifelse(scores$Type == "Del" & scores$`-log10(q-value)` > cutoff, "down", "neutral"))
  # colors to plot
  cnv_palette = c("up"="#bd0000", "down"= "#2e5096", "neutral"="#D2D2D3")
  # if no file is provided, annotate with oncogenes in GAMBLR package
  if(missing(genes_to_label)){
    genes_to_label=GAMBLR::grch37_oncogene %>%
      dplyr::mutate(across(c(chrom, start, end), as.integer)) %>%
      data.table::as.data.table()
  }else{
    genes_to_label=data.table::fread(genes_to_label)
    colnames(genes_to_label)[1:3] <- c("chrom", "start", "end")
    genes_to_label=genes_to_label %>%
      # for now, drop the X chromosome since GISTIC runs without sex chromosmes
      dplyr::filter(!grepl("X", chrom)) %>%
      dplyr::mutate(across(c(chrom, start, end), as.integer)) %>%
      data.table::as.data.table()
  }
  # overlap scores with genes to annotate
  scores <- data.table::foverlaps(scores %>% data.table::setkey(., Chromosome, Start, End),
                      genes_to_label %>% data.table::setkey(., chrom, start, end),
                      by.x=c("Chromosome", "Start", "End"),
                      by.y = c("chrom", "start", "end"),
                      type = "within") %>%
    # if gene to annotate is provided, but it is in region with no CNV, do not label it
    dplyr::mutate(gene=ifelse(!is.na(gene) & fill=="neutral", NA, gene)) %>%
    # if gene is covering multiple adjacent regions, label only once
    dplyr::group_by(gene) %>%
    dplyr::mutate(newcol = ifelse(!is.na(gene) & !duplicated(gene), gene, NA),
           gene = newcol) %>%
    dplyr::select(-newcol)
  # get coordinates to label chromosome numbers
  xses <- scores %>%
    dplyr::group_by(Chromosome) %>%
    dplyr::mutate(End=max(End)/2) %>%
    dplyr::pull(End)

  # main plotting
  ggplot(data = scores,
         aes(x=Start, y=`G-score`,color=fill, label=gene)) +
    geom_bar(size=0.2,stat='identity', position="dodge") +
    ylab('G-score') +
    ggrepel::geom_text_repel(data = subset(scores, !is.na(gene) & Type=="Amp"),
                    nudge_y = max(subset(scores, !is.na(gene) & Type=="Amp")$`G-score`)*adjust_amps,
                    size=label_size,
                    segment.size = 0.5,
                    segment.color = "#000000",
                    force_pull = force_pull,
                    arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
                    segment.curvature = segment.curvature,
                    segment.ncp = segment.ncp,
                    segment.angle = segment.angle) +
    ggrepel::geom_text_repel(data = subset(scores, !is.na(gene) & Type=="Del"),
                             nudge_y = min(subset(scores, !is.na(gene) & Type=="Del")$`G-score`)*adjust_dels,
                             nudge_x=subset(scores, !is.na(gene) & Type=="Del")$Start,
                             size=label_size,
                             segment.size = 0.5,
                             segment.color = "#000000",
                             force_pull = force_pull,
                             arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
                             segment.curvature = segment.curvature,
                             segment.ncp = segment.ncp,
                             segment.angle = segment.angle) +
    facet_grid(. ~ Chromosome, scales="free") +
    scale_color_manual(values=cnv_palette) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=16, colour = "black"),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(colour = "black"),
          legend.position = "none",
          panel.spacing.x=unit(0.1, "lines"),
          panel.border = element_blank(),
          text=element_text(size=16, colour = "black"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.grid = element_blank()) +
    geom_hline(yintercept = 0, size=7) +
    geom_text(aes(label = Chromosome, x = xses, y = 0), size = 4, color="white")
}

#' Define function for consistent plot theme
#'
#' @param base_size Size of the font on the plot. Defaults to 14
#' @param base_family Font family to be used on the plot. Defaults to Arial. Always use cairo device when saving the resulting plot!
#' @param my_legend_position Where to draw the legend? Defaults to the bottom of the plot
#' @param my_legend_direction Which direction to draw the legend? Defaults to horizontal
#'
#'
#' @return nothing
#' @export
#' @import ggplot2 ggthemes
#'
#' @examples
#' ggplot(mpg, aes(displ, hwy, colour = class)) +
#' geom_point() +
#' theme_Morons()

theme_Morons <- function(base_size=14,
                        base_family="Arial",
                        my_legend_position="bottom",
                        my_legend_direction = "horizontal") {
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(colour = "black"),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.2)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = base_size, family=base_family),
            axis.line = element_line(colour="black", size = rel(0.8)),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = my_legend_position,
            legend.direction = my_legend_direction,
            legend.title = element_text(face="italic"),
            strip.background = element_rect(
              color="black", fill="white", size=1, linetype="solid"),
            strip.text = element_text(face="bold")
    ))
}

#' Create a forest plot comparing mutation frequencies for a set of genes between two groups.
#'
#' @param maf A maf data frame. Minimum required columns are Hugo_Symbol and Tumor_Sample_Barcode.
#' @param metadata Metadata for the comparisons. Minimum required columns are Tumor_Sample_Barcode and the column assigning each case to one of two groups.
#' @param comparison_column Mandatory: the name of the metadata column containing the comparison values.
#' @param comparison_values Optional: If the comparison column contains more than two values or is not a factor, specify a character vector of length two in the order you would like the factor levels to be set, reference group first.
#' @param separate_hotspots Optional: If you would like to treat hotspots separately from other mutations in any gene. Requires that the maf file is annotated with GAMBLR::annotate_hotspots.
#' @param comparison_name Optional: Specify the legend title if different from the comparison column name.
#' @param custom_colours Optional: Specify a named vector of colours that match the values in the comparison column.
#' @param custom_labels Optional: Specify custom labels for the legend categories. Must be in the same order as comparison_values.
#' @return A ggplot object with a side-by-side forest plot and bar plot showing mutation incidences across two groups.
#' @export
#' @import dplyr cowplot broom reshape2
#'
#' @examples
#' metadata <- get_gambl_metadata(case_set = "tFL-study") #%>%
#'   dplyr::filter(pairing_status == "matched") %>%
#'   dplyr::filter(consensus_pathology %in% c("FL", "DLBCL"))
#'
#' maf <- get_coding_ssm(limit_samples = metadata$sample_id, basic_columns = TRUE)
#' genes <- c("ATP6V1B2", "EZH2", "TNFRSF14", "RRAGC")
#' comparison_column = "consensus_pathology"
#' comparison_values = c("DLBCL", "FL")
#' comparison_name = "FL vs DLBCL"
#'
#' prettyForestPlot(maf, metadata, genes, comparison_column, comparison_values, separate_hotspots = FALSE, comparison_name)
prettyForestPlot <- function(maf,mutmat, metadata, genes, comparison_column, comparison_values = FALSE, separate_hotspots = FALSE, comparison_name = FALSE, custom_colours = FALSE, custom_labels = FALSE, max_q=1){

  # Subset the maf file to the specified genes
  {
    if(!exists("genes"))
      stop("Please provide a character vector of genes you wish to compare. ")
  }



  # If no comparison_values are specified, derive the comparison_values from the specified comparison_column
  if(comparison_values[1] == FALSE){
    if(class(metadata[[comparison_column]]) == "factor"){
      comparison_values = levels(metadata[[comparison_column]])
    } else {
      comparison_values = unique(metadata[[comparison_column]])
    }
  }

  # Ensure there are only two comparison_values
  {
    if(length(comparison_values) != 2)
      stop(paste0("Your comparison must have two values. \nEither specify comparison_values as a vector of length 2 or subset your metadata so your comparison_column has only two unique values or factor levels."))
  }

  # Subset the metadata to the specified comparison_values and the maf to the remaining sample_ids
  metadata <- metadata[metadata[[comparison_column]] %in% comparison_values, ]


  # Ensure the metadata comparison column is a factor with levels matching the input
  metadata$comparison = factor(metadata[[comparison_column]], levels = comparison_values)

  if(!missing(maf)){
    maf <- maf[maf$Hugo_Symbol %in% genes, ]
    maf <- maf[maf$Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode, ]
  }
  # If separate_hotspots = true, confirm the input maf is hotspot annotated

  if(!missing(mutmat)){
    #add the required columns from the metadata and make the names consistent
    mutmat = left_join(dplyr::select(metadata, sample_id, comparison),mutmat) %>%
      dplyr::rename("Tumor_Sample_Barcode"="sample_id")
  }else if(!missing(maf)){
    if(separate_hotspots){
        if(!"hot_spot" %in% colnames(maf))
          stop("No \"hot_spot\" column in maf file. Annotate your maf file with GAMBLR::annotate_hot_spots() first. ")
      maf$Hugo_Symbol = ifelse(!is.na(maf$hot_spot), paste0(maf$Hugo_Symbol, "_hotspot"), maf$Hugo_Symbol)
    }
    # Convert the maf file to a binary matrix
    mutmat <- maf %>%
      dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      left_join(dplyr::select(metadata, Tumor_Sample_Barcode, comparison),
              by = "Tumor_Sample_Barcode") %>%
      distinct() %>%
      dplyr::mutate(is_mutated = 1) %>%
      pivot_wider(names_from = Hugo_Symbol,
                values_from = is_mutated,
                values_fill = 0) %>%
      dplyr::mutate(across(where(is.numeric), ~replace_na(., 0)))
  }else{
    message("provide a MAF or mutation matrix")
    return()
  }
  fish_test <- mutmat %>%
    pivot_longer(-c(Tumor_Sample_Barcode, comparison),
                 names_to = "gene",
                 values_to = "is_mutated") %>%
    dplyr::mutate(is_mutated = factor(is_mutated, levels = c("1", "0"))) %>%
    group_by(gene) %>%
    dplyr::summarise(table = list(table(is_mutated, comparison))) %>%
    dplyr::mutate(
      test = map(table, fisher.test),
      tidy = map(test, broom::tidy)
    ) %>%
    unnest(tidy) %>%
    dplyr::mutate(q.value = p.adjust(p.value, "BH")) %>%
    dplyr::select(-c(table, test, method, alternative)) %>%
    dplyr::filter(q.value <= max_q) %>%
    dplyr::mutate(gene = fct_reorder(gene, estimate))
  #fish_test <- mutate(fish_test, gene = fct_reorder(gene, estimate))
  #fish_test$gene = factor(fish_test$gene,levels=unique(fish_test$gene))

  point_size = 50/round(length(fish_test$gene))
  if(point_size<1){
    point_size = 1
  }
  font_size = 360/round(length(fish_test$gene))
  if(font_size<4){
    font_size=4
  }else if(font_size > 20){
    font_size = 20
  }
  message(paste("FONT:",font_size,"POINT:",point_size,length(fish_test$gene)))
  forest <- fish_test %>%
    ggplot(aes(x = gene, y = log(estimate))) +
    geom_point(size = point_size, shape = "square") +
    geom_hline(yintercept = 0, lty = 2) +
    coord_flip() +
    geom_errorbar(aes(ymin = log(conf.low), ymax = log(conf.high), width = 0.2)) +
    ylab("ln(Odds Ratio)") +
    xlab("Mutated Genes") +
    cowplot::theme_cowplot() +
    theme(axis.text.y = element_text(size = font_size))

  if(comparison_name == FALSE){
    comparison_name = comparison_column
  }

  if(custom_colours[1] == FALSE){
    if(length(levels(metadata$comparison)[levels(metadata$comparison) %in% names(get_gambl_colours())]) == 2){
      colours = get_gambl_colours()[levels(metadata$comparison)]
    } else {
      colours = get_gambl_colours(classification = "blood")[c("Red", "Blue")]
      names(colours) = levels(metadata$comparison)
    }
  } else {
    colours <- custom_colours
  }

  if(custom_labels[1] == FALSE){
    labels = levels(metadata$comparison)
    names(labels) = levels(metadata$comparison)
  } else if(length(custom_labels) != 2) {
    labels = levels(metadata$comparison)
    names(labels) = levels(metadata$comparison)
    print("Provided custom labels is not a character vector of length 2. Defaulting to comparison factor levels as labels. ")
  } else {
    labels = custom_labels
    names(labels) = comparison_values
  }

  bar <- mutmat %>%
    dplyr::select(-Tumor_Sample_Barcode) %>%
    reshape2::melt(., id.vars = c("comparison"), value.name="is_mutated", variable.name="gene") %>%
    group_by(gene, comparison) %>%
    drop_na() %>%
    summarise(percent_mutated = sum(is_mutated)/n() * 100) %>%
    dplyr::filter(gene %in% fish_test$gene) %>%
    dplyr::mutate(gene = factor(gene, levels = levels(fish_test$gene))) %>%
    ggplot(aes(x = gene, y = percent_mutated, fill = comparison)) +
    geom_col(position = "dodge", width = 0.5) +
    xlab("") + ylab("% Mutated") +
    coord_flip() +
    scale_fill_manual(name = comparison_name, values = colours, labels = labels[levels(metadata$comparison)]) +
    cowplot::theme_cowplot() +
    theme(axis.text.y = element_blank(),
          legend.position = "bottom",
          legend.justification = )

  legend = cowplot::get_legend(bar)

  plots <- plot_grid(forest, bar + theme(legend.position = "none"), rel_widths = c(1, 0.6), nrow = 1)

  arranged_plot = cowplot::plot_grid(plot_grid(NULL, legend, NULL, nrow = 1), plots, nrow = 2, rel_heights = c(0.1, 1))

  return(list(fisher=fish_test,forest=forest,bar=bar,legend=legend,arranged=arranged_plot))
}

#' Make an heatmap that is looking cute using ComplexHeatmap. The metadata is expected to follow the structure and column naming used in GAMBL.
#' If you provide your own non-GAMBL samples and metadata, you must include at least the columns with names corresponding to annotation tracks and column "Tumor_Sample_Barcode"
#' showing sample ids. The metadata can contain numeric columns, which will be plotted as numeric variables in the annotation. The efature matrix is supplied in this_matrix argument
#' and is expected to have samples in rows, and features in columns. The argument importance_values is similar to the widths of NMF object or importance values for feature/group from RF models.
#' It is also expected to have column names (having names of the groups that will be shown on heatmap) and rownames (corresponding to feature ids).
#' @param this_matrix A data frame with column Tumor_Sample_Barcode and a column for each feature. Can be binary. Expected to not contain negative values.
#' @param importance_values Provide a data frame of feature (in rows) by group (in columns) with numeric values representative of feature importance. Can be obtained from rf$inportance or basis(NMF)
#' @param these_samples_metadata Data frame containing metadata for your samples
#' @param max_number_of_features_per_group Optional argument to indicate how many features from each group to be considered for display. Default is 10
#' @param splitColumnName Optional argument to indicate which metadata column to split on. Default is set to pathology
#' @param metadataColumns A vector containing the categorical column names you want to plot below
#' @param numericMetadataColumns A vector containing the numeric columns you want to plot below
#' @param numericMetadataMax A numeric vector of cutoffs to apply to numeric columns above
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette
#' @param legend_direction Optional argument to indicate whether legend should be in horizontal (default) or vertical position
#' @param legend_position Optional argument to indicate where the legend should be drawn. The default is set to bottom, but can also accept top, right, and left.
#' @param legend_row Fiddle with these to widen or narrow your legend (default 3)
#' @param legend_col Fiddle with these to widen or narrow your legend (default 3)
#' @param fontSizeGene Font size for gene labels (default 6)
#' @param metadataBarHeight Optional argument to adjust the height of bar with annotations. The default is 1.5
#' @param leftStackedWidth Optional argument to control how wide should the stacked plot on the left be. The default is 4
#' @param metadataBarFontsize Optional argument to control for the font size of metadata annotations. The default is 5
#' @param groupNames optional vector of group names to be displayed above heatmap. Should be the same length as the number of groups that will be shown. Default is NULL (no labels)
#'
#' @return
#' @export
#' @import ComplexHeatmap grid dplyr circlize
#'
#' @examples
#' splendidHeatmap(
#'  this_matrix = data,
#'  importance_values = rf$importance[,c(1:3)],
#'  these_samples_metadata = MASTER.METADATA,
#'  splitColumnName = "pathology",
#'  metadataColumns = c("cohort", "pathology", "sex", ".", "COO_consensus", "DHITsig_consensus", "seq_type"),
#'  numericMetadataColumns = ".",
#'  numericMetadataMax = 0.7,
#'  custom_colours=custom_colours)

splendidHeatmap = function(this_matrix,
                           importance_values,
                           these_samples_metadata,
                           max_number_of_features_per_group = 10,
                           splitColumnName = "pathology",
                           metadataColumns = c("pathology"),
                           numericMetadataColumns = NULL,
                           numericMetadataMax = NULL,
                           custom_colours=NULL,
                           legend_direction="horizontal",
                           legend_position="bottom",
                           legend_row=3,
                           legend_col=3,
                           fontSizeGene=6,
                           metadataBarHeight=1.5,
                           leftStackedWidth=4,
                           metadataBarFontsize=5,
                           groupNames = NULL){
  
  comparison_groups <- unique(these_samples_metadata[,splitColumnName])

  if(!is.null(splitColumnName) & (splitColumnName %in% metadataColumns)){
    metadataColumns <- c(splitColumnName, metadataColumns[!metadataColumns==splitColumnName])
  }

  if(!is.null(numericMetadataColumns) & length(intersect(numericMetadataColumns, metadataColumns))>0){
    message(paste0("The column(s) ", numericMetadataColumns, " specified both in metadata and numeric metadata. Plotting as numeric values..."))
    metadataColumns = metadataColumns[!metadataColumns %in% numericMetadataColumns]
  }

  # get which group samples belong to
  metadata_df <- these_samples_metadata[,c("Tumor_Sample_Barcode", metadataColumns, numericMetadataColumns)] %>%
    as.data.frame() %>%
    column_to_rownames(., "Tumor_Sample_Barcode")

  if(!is.null(numericMetadataMax)){
      max_list <- setNames(numericMetadataMax,numericMetadataColumns)
      metadata_df <- metadata_df %>%
        dplyr::mutate(across(names(max_list), ~ ifelse(.x > max_list[[cur_column()]], max_list[[cur_column()]], .x)))
  }

  my_colours <- NULL
  these_names=NULL
  for (i in 1:length(metadataColumns)){
    this_metadata_column <- get_gambl_colours(metadataColumns[i])
    if (sum(is.na(names(this_metadata_column[unique(these_samples_metadata[,metadataColumns[i]])])))<=1 &
        length(unique(these_samples_metadata[,metadataColumns[i]])) > 1){
      these_names = c(these_names,metadataColumns[i])
      my_colours = append(my_colours, list(c(this_metadata_column, "NA"="#BDBDC1FF")))
      names(my_colours) = these_names
    }
  }

  my_colours <- c(custom_colours, my_colours)

  col_fun=circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  for(exp in numericMetadataColumns){
    my_colours[[exp]] = col_fun
  }

  # get all features
  w <- importance_values[,comparison_groups]
  w <- as.data.frame(w) %>%
    dplyr::mutate_if(is.character,as.numeric)



  # extract most important features, while taking the feature with highest weight for a particular cluster if it was seen before for other cluster with lower weight
  FEATURES <- w[,1] %>%
    as.data.frame() %>% 
    `rownames<-`(rownames(w)) %>%
    dplyr::arrange(desc(.)) %>%
    head(., max_number_of_features_per_group) %>%
    rownames_to_column(., var="Feature") %>%
    dplyr::mutate(group=comparison_groups[1])
  for (i in 2:length(comparison_groups)){
    FEATURES <- rbind(as.data.frame(FEATURES),
                    w[,i] %>% as.data.frame() %>%
                      `rownames<-`(rownames(w)) %>%
                      dplyr::arrange(desc(.)) %>%
                      head(., max_number_of_features_per_group+3) %>%
                      rownames_to_column(., var="Feature") %>%
                      dplyr::mutate(group=comparison_groups[i])) %>%
    dplyr::group_by(Feature) %>%
    dplyr::filter(. == max(.)) %>%
    dplyr::arrange(group)
  }
  FEATURES <- as.data.frame(FEATURES)

  mat <- this_matrix %>%
    merge(., metadata_df %>%
             rownames_to_column(., "Tumor_Sample_Barcode") %>%
             dplyr::select(Tumor_Sample_Barcode, splitColumnName)) %>%
             as.data.frame()
  mat[,splitColumnName] = factor(mat[,splitColumnName])

  # breaks used to display groups with different colors on heatmap
  bk <- c(0,seq(0.5, length(comparison_groups)+0.5, 1))

  # colors used to show on the heatmap body. Starts with white - the color of feature absence
  my_palette <- c("white", rev(unlist(my_colours[splitColumnName])))
  my_palette <- unname(my_palette)

  # get each group and label the events for each feature with group number
  mat_2 <- mat[,-ncol(mat)]
  # subset samples of each group
  MY.LIST <- list()
  for (i in 1:length(comparison_groups)){
    MY.LIST[[i]] <- assign(comparison_groups[i], mat_2 %>%
                           as.data.frame(.) %>%
                           column_to_rownames(., var="Tumor_Sample_Barcode") %>%
                           t(.) %>%
                           as.data.frame(.) %>%
                           dplyr::select(metadata_df %>%
                                    dplyr::filter(base::get(splitColumnName)==comparison_groups[i]) %>%
                                    rownames) )
  }

  # assign numbers - used for coloring of heatmap body
  for(i in 1:length(comparison_groups)){
    MY.LIST[[i]][MY.LIST[[i]]>0] <- i
  }

  # bind them all together for plotting
  mat_2 <- do.call(cbind, MY.LIST) %>%
    as.data.frame(.) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(., var="Tumor_Sample_Barcode") %>%
    base::merge(., metadata_df %>%
                rownames_to_column(., "Tumor_Sample_Barcode") %>%
                dplyr::select(Tumor_Sample_Barcode, splitColumnName),
                by="Tumor_Sample_Barcode")


  # specify where row breaks should be on heatmap
  breaks <- 0
  for (this_group in comparison_groups){
    N <- (nrow(FEATURES %>% dplyr::filter(group==this_group)))
    breaks <- c(breaks, N)
  }

  # second, make a vector that will be supplied to ComplexHeatmap
  my_vector <- NULL
  for (i in 1:(length(breaks))){
    my_vector <- c(my_vector,
                 rep(i-1, breaks[i]))
  }

  # prepare matrix for stacked barplots on the left
  STACKED <- data.frame(matrix(NA, ncol=1, nrow=nrow(FEATURES)))[-1]
  rownames(STACKED) <- FEATURES$Feature
  for (i in 1:length(comparison_groups)) {
  STACKED <- cbind(STACKED,
                   mat_2[,c("Tumor_Sample_Barcode", FEATURES$Feature)] %>%
                     base::merge(., metadata_df %>%
                                   rownames_to_column(., "Tumor_Sample_Barcode") %>%
                                   dplyr::select(Tumor_Sample_Barcode, splitColumnName),
                                   by="Tumor_Sample_Barcode") %>%
                     dplyr::arrange(!!sym(splitColumnName)) %>%
                     dplyr::filter(base::get(splitColumnName)==comparison_groups[i]) %>%
                     dplyr::select(-Tumor_Sample_Barcode, -splitColumnName) %>%
                     dplyr::summarise_all(funs(sum)) %>%
                     t(.) %>%
                     `colnames<-`(comparison_groups[i]) %>%
                     as.data.frame(.) %>%
                     dplyr::mutate_all(~(./i)/nrow(metadata_df)))
  }
  m <- t(apply(STACKED, 1, function(x) x/sum(x)))

  used_for_ordering_df <- t(base::merge(mat_2 %>%
                                            dplyr::select(-splitColumnName),
                                        metadata_df %>%
                                            rownames_to_column(., "Tumor_Sample_Barcode"),
                                        by="Tumor_Sample_Barcode") %>%
                          column_to_rownames(., var="Tumor_Sample_Barcode") %>%
                          dplyr::arrange(!!!syms(metadataColumns), desc(!!!syms(numericMetadataColumns))) %>%
    dplyr::select(FEATURES$Feature))
  
  used_for_ordering <- colnames(used_for_ordering_df)

  # left annotation: stacked feature weights
  ha = rowAnnotation(`feature abundance` = anno_barplot(m, gp = gpar(fill = my_palette[1:length(comparison_groups)+1]),
                                                      bar_width = 1, width = unit(leftStackedWidth, "cm"), 
                                                      axis_param = list(side = legend_position, labels_rot = 0)))

  # bottom annotation: tracks indicating metadata
  ha_bottom = HeatmapAnnotation(df = metadata_df[ (order(match(rownames(metadata_df), used_for_ordering))), ] %>%
                                dplyr::arrange(!!!syms(metadataColumns), desc(!!!syms(numericMetadataColumns))) %>%
                                dplyr::select(-splitColumnName),
                              col = my_colours,
                              simple_anno_size = unit(metadataBarHeight, "mm"),
                              gap = unit(0.25*metadataBarHeight, "mm"),
                              annotation_name_gp=gpar(fontsize=metadataBarFontsize),
                              annotation_legend_param =
                                list(nrow=legend_row,
                                     ncol=legend_col,
                                     direction=legend_direction))

  # top annotation: groups of interest to split on
  ha_top = HeatmapAnnotation(df = metadata_df[ (order(match(rownames(metadata_df), used_for_ordering))), ] %>%
                             dplyr::arrange(!!!syms(metadataColumns), desc(!!!syms(numericMetadataColumns))) %>%
                             dplyr::select(splitColumnName),
                           col = my_colours[splitColumnName],
                           simple_anno_size = unit(metadataBarHeight, "mm"),
                           gap = unit(0.25*metadataBarHeight, "mm"),
                           annotation_name_gp=gpar(fontsize=fontSizeGene*1.5),
                           annotation_legend_param =
                             list(nrow=legend_row,
                                  ncol=legend_col,
                                  direction=legend_direction))

  splendidHM <- ComplexHeatmap::Heatmap(used_for_ordering_df,
                               col = my_palette,
                               show_column_names = FALSE,
                               cluster_columns = FALSE,
                               cluster_rows = FALSE,
                               row_names_gp = gpar(fontsize = fontSizeGene),
                               show_heatmap_legend = FALSE,
                               row_split = my_vector,
                               row_title = NULL,
                               left_annotation=ha,
                               bottom_annotation=ha_bottom,
                               top_annotation=ha_top,
                               column_split=dplyr::pull(metadata_df[(order(match(rownames(metadata_df), used_for_ordering))), ], splitColumnName),
                               column_title=groupNames)

  draw(splendidHM, heatmap_legend_side = legend_position, annotation_legend_side = legend_position)

}
